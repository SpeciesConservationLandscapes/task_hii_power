import argparse
import ee
import os
from task_base import HIITask


CALC_PREVIOUS_ANNUAL_VIIRS = "viirs"
CALC_CALIBRATION_COEFFICIENTS = "coefficients"
CALC_WEIGHTING_QUANTILES = "quantiles"
JOB_CHOICES = (
    CALC_PREVIOUS_ANNUAL_VIIRS,
    CALC_CALIBRATION_COEFFICIENTS,
    CALC_WEIGHTING_QUANTILES,
)


class HIIPowerPreprocessTask(HIITask):
    inputs = {
        "dmsp_viirs_calibrated": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/source/nightlights/dmsp_viirs_calibrated",
            "static": True,
        },
        "viirs": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG",
            "maxage": 1,
        },
        "watermask": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True,
        },
    }

    thresholds = {
        "latitude": {"min_lat": 0, "min_val": 0.1, "max_lat": 60, "max_val": 0.75},
    }
    calibration_coefficients = {"slope": 10.53, "intercept": 24.62}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.job = (
            kwargs.get("job") or os.environ.get("job") or CALC_PREVIOUS_ANNUAL_VIIRS
        )

        self.dmsp = ee.ImageCollection(
            self.inputs["dmsp_viirs_calibrated"]["ee_path"]
        ).filterDate("1992", "2013")
        self.viirs = ee.ImageCollection(self.inputs["viirs"]["ee_path"])
        self.watermask = ee.Image(self.inputs["watermask"]["ee_path"])

    def viirs_annual_image(self, year):
        year_int = ee.Number(year).int()
        year_str = year_int.format("%d")
        start_date = ee.Date.parse("YYYY", year_str)
        end_date = start_date.advance(1, "year")

        min_lat = self.thresholds["latitude"]["min_lat"]
        max_lat = self.thresholds["latitude"]["max_lat"]
        min_val = self.thresholds["latitude"]["min_val"]
        max_val = self.thresholds["latitude"]["max_val"]

        def clean_images(image):
            viirs_radiance = image.select("avg_rad")
            cleaned = viirs_radiance.updateMask(viirs_radiance.gt(0))
            return cleaned

        viirs_annual_median = (
            self.viirs.filterDate(start_date, end_date).map(clean_images).median()
        )

        latitude = (
            ee.Image.pixelLonLat()
            .select("latitude")
            .abs()
            .pow(4)
            .clamp(min_lat, max_lat ** 4)
        )

        thresholds = (
            (latitude.divide(max_lat ** 4).multiply(max_val - min_val).add(min_val))
            .clamp(min_val, max_val)
            .rename("latitude_thresholds")
        )

        annual_viirs_mask = viirs_annual_median.multiply(0.25)
        annual_viirs = viirs_annual_median.mask(annual_viirs_mask)
        latitude_threshold_mask = annual_viirs.mask().gte(thresholds)
        annual_viirs = (
            annual_viirs.updateMask(latitude_threshold_mask)
            .rename("viirs")
            .set({"system:time_start": start_date.millis()})
        )
        return annual_viirs

    def stable_lights_image(self, std_dev):
        dmsp_recent = self.dmsp.filterDate("2000", "2013")
        dmsp_std_dev = dmsp_recent.reduce(ee.Reducer.stdDev())
        dmsp_stable = dmsp_std_dev.lte(std_dev)
        return dmsp_stable

    def stable_light_points(self):
        stable_lights = self.stable_lights_image(2)
        dmsp_latest = self.dmsp.filterDate("2012", "2013").first()
        dmsp_stable = (
            dmsp_latest.updateMask(stable_lights)
            .updateMask(self.watermask)
            .updateMask(dmsp_latest.gt(0))
            .rename("DMSP")
        )
        viirs_earliest = (
            self.viirs_annual_image(2014).updateMask(dmsp_stable).rename("VIIRS")
        )

        viirs_earliest = viirs_earliest.updateMask(viirs_earliest.gt(0))

        dmsp_class = dmsp_stable.round().int().rename("NL_CLASS")
        regression_image = dmsp_class.addBands([dmsp_stable, viirs_earliest])
        region = ee.Geometry.Polygon(self.extent, proj=self.crs, geodesic=False)

        return regression_image.stratifiedSample(
            numPoints=5,
            classBand="NL_CLASS",
            region=region,
            scale=self.scale,
            projection=self.crs,
            tileScale=16,
            dropNulls=True,
        )

    def viirs_to_dmsp(self):
        viirs_projection = self.viirs.first().projection()
        dmsp_projection = self.dmsp.first().projection()
        annual_viirs = (
            self.viirs_annual_image(self.taskdate.year)
            .log()
            .multiply(self.calibration_coefficients["slope"])
            .add(self.calibration_coefficients["intercept"])
            .unmask(0)
            .setDefaultProjection(viirs_projection)
            .reduceResolution(ee.Reducer.mean())
            .reproject(dmsp_projection)
            .updateMask(self.watermask)
            .clamp(0, 63)
            .uint8()
            .selfMask()
        )
        return annual_viirs

    def quantile_calc(self, image_collection, year):
        year_int = ee.Number(year).int()
        year_str = year_int.format("%d")
        start_date = ee.Date.parse("YYYY", year_str)
        end_date = start_date.advance(1, "year")

        nightlights_year = image_collection.filterDate(start_date, end_date).first()
        region = ee.Geometry.Polygon(self.extent, proj=self.crs, geodesic=False)
        projection = image_collection.projection()
        scale = projection.nominalScale()
        nightlight_quantiles = nightlights_year.reduceRegion(
            reducer=ee.Reducer.percentile([10, 20, 30, 40, 50, 60, 70, 80, 90]),
            geometry=region,
            crs=projection,
            scale=scale,
            maxPixels=self.ee_max_pixels,
            tileScale=16,
        )

        return nightlight_quantiles

    def calc(self):
        if self.job == CALC_PREVIOUS_ANNUAL_VIIRS:
            # Running this task for any date in (e.g.) 2020 -- even 2020-01-01 -- will produce
            # an annual calibrated_viirs from 2020-01-01 through 2020-12-31
            exportpath = "source/nightlights/dmsp_viirs_calibrated"
            calibrated_viirs = self.viirs_to_dmsp()
            self.export_image_ee(calibrated_viirs, exportpath, self.taskdate.year)

        elif self.job == CALC_CALIBRATION_COEFFICIENTS:
            regression_points = self.stable_light_points()
            # TODO: operationalize in pandas this R code for calculating regression coefficients
            """
            regression_points = regression_points[
                which(regression_points$VIIRS>0 &
                regression_points$DMSP>0),
            ]

            regression = lm(regression_points$DMSP ~ log(regression_points$VIIRS))

            slope = regression$coefficients[1]
            intercept = regression$coefficients[1]
            r2 = summary(regression)[8]
            r2_adj = summary(regression)[9]
            print('slope: ', slope, '\n intercept: ', intercept, '\n r2: ', r2, '\n r2_adj: ', r2_adj)
            """

        elif self.job == CALC_WEIGHTING_QUANTILES:
            quantiles = self.quantile_calc(self.dmsp, 1992)
            print(quantiles.getInfo())

    def check_inputs(self):
        super().check_inputs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate")
    # TODO: in README, document usage of this task with custom args separately and called from main task
    parser.add_argument(
        "-j",
        "--job",
        nargs="?",
        choices=JOB_CHOICES,
        default=CALC_PREVIOUS_ANNUAL_VIIRS,
        const=CALC_PREVIOUS_ANNUAL_VIIRS,
        help="Calculation to run with HIIPowerPreprocessTask.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing outputs instead of incrementing",
    )
    options = parser.parse_args()
    power_task = HIIPowerPreprocessTask(**vars(options))
    power_task.run()
