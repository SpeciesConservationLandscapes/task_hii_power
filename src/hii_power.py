import argparse
import ee
from datetime import datetime, timezone
from task_base import EETask


class HIIPower(EETask):
    ee_rootdir = "projects/HII/v1/sumatra_poc"
    ee_driverdir = "driver/power"
    inputs = {
        "viirs": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": "NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG",
        },
        "watermask": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
        },
    }
    scale = 300
    gpw_cadence = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))

    def calc(self):
        viirs = ee.ImageCollection(self.inputs["viirs"]["ee_path"])
        watermask = ee.Image(self.inputs["watermask"]["ee_path"])
        polar = ee.Feature(
            ee.FeatureCollection(
                "projects/HII/v1/misc/polar_for_viirs_correction"
            ).first()
        )

        # TODO: add property to source fc to avoid this step (and "foo"). Also, better to use server-side List()
        # instead of python list as argument to FeatureCollection()?
        viirs_polar_correction = ee.FeatureCollection(
            [polar.set({"foo": 1})]
        ).reduceToImage(["foo"], ee.Reducer.sum())

        def viirs_func(year):
            startdate = ee.Date.fromYMD(year, 1, 1)
            viirs_median = (
                viirs.filterDate(startdate, startdate.advance(1, "year"))
                .select("avg_rad")
                .median()
                .reproject(crs=self.crs, scale=self.scale)
            )

            viirs_outside_solar = viirs_median.multiply(viirs_polar_correction.neq(1))
            viirs_inside_solar = viirs_median.multiply(viirs_polar_correction.eq(1))
            viirs_no_solar = viirs_outside_solar.add(viirs_inside_solar)
            return viirs_no_solar

        # TODO: replace hardcoded arguments with class properties
        def calibration(image):
            image = image.unmask()
            return ee.Image(1569.44).add(ee.Image(828.19).multiply(image.log()))

        def reproject_dmsp(image):
            return image.reproject(crs=self.crs, scale=self.scale)

        # TODO: These should be moved into an image collection, added to self.inputs, and attributed
        # with system:time_start and a maxage property = 1.
        # Then merge this dmsp ic with the viirs ic, refactor code to use self.get_most_recent_image()
        # So most recent year (DMSP) or most recent month (VIIRS), and then we don't have to
        # run viirs_func over all viirs images, and in fact can avoid that median step altogether
        # (unless we want to use the median from one time point to either side?)
        dmsp_1992 = ee.Image("users/aduncan/yale/F101992")
        dmsp_1993 = ee.Image("users/aduncan/yale/F101993")
        dmsp_1994 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F101994"),
                ee.Image("users/aduncan/yale/F121994"),
            ]
        ).mean()
        dmsp_1995 = ee.Image("users/aduncan/yale/F121995")
        dmsp_1996 = ee.Image("users/aduncan/yale/F121996")
        dmsp_1997 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F121997"),
                ee.Image("users/aduncan/yale/F141997"),
            ]
        ).mean()
        dmsp_1998 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F121998"),
                ee.Image("users/aduncan/yale/F141998"),
            ]
        ).mean()
        dmsp_1999 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F121999"),
                ee.Image("users/aduncan/yale/F141999"),
            ]
        ).mean()
        dmsp_2000 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F142000"),
                ee.Image("users/aduncan/yale/F152000"),
            ]
        ).mean()
        dmsp_2001 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F142001"),
                ee.Image("users/aduncan/yale/F152001"),
            ]
        ).mean()
        dmsp_2002 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F142002"),
                ee.Image("users/aduncan/yale/F152002"),
            ]
        ).mean()
        dmsp_2003 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F142003"),
                ee.Image("users/aduncan/yale/F152003"),
            ]
        ).mean()
        dmsp_2004 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F152004"),
                ee.Image("users/aduncan/yale/F162004"),
            ]
        ).mean()
        dmsp_2005 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F152005"),
                ee.Image("users/aduncan/yale/F162005"),
            ]
        ).mean()
        dmsp_2006 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F152006"),
                ee.Image("users/aduncan/yale/F162006"),
            ]
        ).mean()
        dmsp_2007 = ee.ImageCollection(
            [
                ee.Image("users/aduncan/yale/F152007"),
                ee.Image("users/aduncan/yale/F162007"),
            ]
        ).mean()
        dmsp_2008 = ee.Image("users/aduncan/yale/F162008")
        dmsp_2009 = ee.Image("users/aduncan/yale/F162009")
        dmsp_2010 = ee.Image("users/aduncan/yale/F182010")
        dmsp_2011 = ee.Image("users/aduncan/yale/F182011")
        dmsp_2012 = ee.Image("users/aduncan/yale/F182012")

        # TODO: see above. This shouldn't be future-constrained.
        yearlist = ee.List.sequence(2014, 2019)

        """
        viirs_ic = ee.ImageCollection(yearlist.map(function(year){return viirs_func(year)}))
          .map(function(image){return calibration(image)});
        """
        viirs_ic = ee.ImageCollection(yearlist.map(viirs_func)).map(calibration)

        # TODO: interpolate 2013 between DMSP 2012 and VIIRS 2014
        dmsp_ic = ee.ImageCollection(
            [
                dmsp_1992,
                dmsp_1993,
                dmsp_1994,
                dmsp_1995,
                dmsp_1996,
                dmsp_1997,
                dmsp_1998,
                dmsp_1999,
                dmsp_2000,
                dmsp_2001,
                dmsp_2002,
                dmsp_2003,
                dmsp_2004,
                dmsp_2005,
                dmsp_2006,
                dmsp_2007,
                dmsp_2008,
                dmsp_2009,
                dmsp_2010,
                dmsp_2011,
                dmsp_2012,
                dmsp_2012,
            ]
        ).map(reproject_dmsp)

        # NEED TO SORT DESCENDING!!! This was a map but don't need right now...
        ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        # TODO: remove hardcoded 1992, use class methods for date validation
        year_num = ee.Number(ee_taskdate.get("year")).subtract(1992)

        # TODO: see above on get_most_recent_image
        hii_power_oneyear = ee.Image(dmsp_ic.merge(viirs_ic).toList(30).get(year_num))

        # TODO: replace hardcoded arguments with class properties
        # TODO: could combine two multiply statements -- but depends on what they mean/do
        hii_power_driver = (
            hii_power_oneyear.unmask(0)
            .multiply(0.01)
            .updateMask(watermask)
            .multiply(10)
        )
        self.export_image_ee(
            hii_power_driver, "{}/{}".format(self.ee_driverdir, "hii_power_driver")
        )

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--taskdate", default=datetime.now(timezone.utc).date())
    options = parser.parse_args()
    power_task = HIIPower(**vars(options))
    power_task.run()
