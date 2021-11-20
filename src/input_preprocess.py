import argparse
import ee
from datetime import datetime
from task_base import HIITask


class HII_Power_Input_Processing(HIITask):
    scale = 300

    inputs = {
        "dmsp_viirs_calibrated": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/source/nightlights/dmsp_viirs_calibrated",
            "maxage": 1,
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
        "latitude": {"min_lat": 0, "min_val": 0.1, "max_lat": 60, "max_val": 0.75,},
    }
    calibration_coefficients = {
        "slope": 10.53,
        "intercept": 24.62,
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dmsp = ee.ImageCollection(
            self.inputs["dmsp_viirs_calibrated"]["ee_path"]
        ).filterDate("1992", "2013")

        self.viirs = ee.ImageCollection(self.inputs["viirs"]["ee_path"])

        self.watermask = ee.Image(self.inputs["watermask"]["ee_path"])

    def viirs_annual_image(
        self, year, viirs_collection, min_lat, max_lat, min_val, max_val
    ):

        year_int = ee.Number(year).int()
        year_str = year_int.format("%d")

        start_date = ee.Date.parse("YYYY", year_str)
        end_date = start_date.advance(1, "year")

        def clean_images(image):
            viirs_radiance = image.select("avg_rad")
            cleaned = viirs_radiance.updateMask(viirs_radiance.gt(0))
            return cleaned

        viirs_annual_median = (
            viirs_collection.filterDate(start_date, end_date).map(clean_images).median()
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

    def stable_lights_image(self, dmsp_collection, std_dev):
        dmsp_recent = dmsp_collection.filterDate("2000", "2013")
        dmsp_std_dev = dmsp_recent.reduce(ee.Reducer.stdDev())
        dmsp_stable = dmsp_std_dev.lte(std_dev)
        return dmsp_stable

    def stable_light_points(
        self, dmsp_collection, viirs_collection, stable_lights_image, watermask
    ):
        dmsp_latest = dmsp_collection.filterDate("2012", "2013").first()
        dmsp_stable = (
            dmsp_latest.updateMask(stable_lights_image)
            .updateMask(watermask)
            .updateMask(dmsp_latest.gt(0))
            .rename("DMSP")
        )
        viirs_earliest = (
            self.viirs_annual_image(
                2014,
                viirs_collection,
                self.thresholds["latitude"]["min_lat"],
                self.thresholds["latitude"]["max_lat"],
                self.thresholds["latitude"]["min_val"],
                self.thresholds["latitude"]["max_val"],
            )
            .updateMask(dmsp_stable)
            .rename("VIIRS")
        )

        viirs_earliest = viirs_earliest.updateMask(viirs_earliest.gt(0))

        dmsp_class = dmsp_stable.round().int().rename("NL_CLASS")
        regression_image = dmsp_class.addBands([dmsp_stable, viirs_earliest])

        # TODO: global bounds?
        gBounds = ee.Geometry.Polygon(
            [-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], None, False
        )

        regression_sample = regression_image.stratifiedSample(
            numPoints=5,
            classBand="NL_CLASS",
            region=gBounds,
            scale=self.scale,
            projection=self.crs,
            tileScale=16,
            dropNulls=True,
        )

    def viirs_to_dmsp(self, viirs_collection, year, viirs_proj, dmsp_proj, watermask):
        annual_viirs = (
            self.viirs_annual_image(
                year,
                viirs_collection,
                self.thresholds["latitude"]["min_lat"],
                self.thresholds["latitude"]["max_lat"],
                self.thresholds["latitude"]["min_val"],
                self.thresholds["latitude"]["max_val"],
            )
            .log()
            .multiply(self.calibration_coefficients["slope"])
            .add(self.calibration_coefficients["intercept"])
            .unmask(0)
            .setDefaultProjection(viirs_proj)
            .reduceResolution(ee.Reducer.mean())
            .reproject(dmsp_proj)
            .updateMask(watermask)
            .clamp(0, 63)
            .uint8()
            .selfMask()
        )
        return annual_viirs

    def calc(self):
        viirs_projection = self.viirs.first().projection()
        dmsp_projection = self.dmsp.first().projection()

        stable_lights = self.stable_lights_image(self.dmsp, 2)
        regression_points = self.stable_light_points(
            self.dmsp, self.viirs, stable_lights, self.watermask
        )
        # TODO: export the regression sample for analysis
        # export to GCS as CSV and then bring back in to run regression in this task?

        calibrated_viirs = self.viirs_to_dmsp(
            self.viirs,
            self.taskdate.year,
            viirs_projection,
            dmsp_projection,
            self.watermask,
        )

        #   var exportString = exportYear.toString();
        #   var id = 'projects/HII/v1/source/nightlights/dmsp_viirs_calibrated/Harmonized_DN_NTL_' + exportString + '_calDMSP';
        #
        #   var description = 'Harmonized_DN_NTL_' + exportString + '_calDMSP';
        #
        #   var dmsp_latest = li_et_al.filterDate('2012', '2013').first()
        #     // .updateMask(watermask)
        #     .rename('DMSP');
        #   // Map.addLayer(dmsp_latest)
        #   var dmsp_1992 = li_et_al.filterDate('1992', '1993').first()
        #     .updateMask(watermask)
        #     .selfMask()
        #     .rename('DMSP');
        #   dmsp_1992 = dmsp_1992.updateMask(dmsp_1992.gt(1));
        #
        #   var viirs_callibrated =  annual_viirs(viirs, year, false)//.select(0)
        #   if (exportYear == 2013) {
        #     var date = ee.Date.parse('YYYY', '2013').millis()
        #     viirs_callibrated = viirs_callibrated.add(dmsp_latest)
        #       .divide(2)
        #       .round()
        #       .uint8()
        #       .set({'system:time_start': date})
        #   } else {
        #     viirs_callibrated = viirs_callibrated
        #   }
        # }

        # self.export_image_ee(hii_power_driver, "driver/power")

    def check_inputs(self):
        super().check_inputs()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing outputs instead of incrementing",
    )
    options = parser.parse_args()
    power_task = HII_Power_Input_Processing(**vars(options))
    power_task.run()
