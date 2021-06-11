import argparse
import ee
from datetime import datetime
from task_base import HIITask


class HIIPower(HIITask):
    DMSP_ERA = datetime(2013, 1, 1)
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
            "maxage": 2,
        },
        "watermask": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True,
        },
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dmsp_viirs_calibrated = ee.ImageCollection(
            self.inputs["dmsp_viirs_calibrated"]["ee_path"]
        )
        self.viirs = ee.ImageCollection(self.inputs["viirs"]["ee_path"])
        self.watermask = ee.Image(self.inputs["watermask"]["ee_path"])
        self.nightlights = None
        self.nightlights_date = None

    def populate_nightlights(self):
        ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        self.nightlights, self.nightlights_date = self.get_most_recent_image(
            self.dmsp_viirs_calibrated
        )
        age = ee_taskdate.difference(self.nightlights_date, "year").getInfo()
        if age > self.inputs["dmsp_viirs_calibrated"]["maxage"]:
            # export viirs images here and set self.nightlights, self.nightlights_date
            self.wait()

    def calc(self):
        self.populate_nightlights()
        # use self.nightlights to create driver

        # hii_power_driver = nightlights.multiply(0.01).updateMask(self.watermask)
        # self.export_image_ee(
        #     hii_power_driver, f"driver/power"
        # )

    def check_inputs(self):
        if self.taskdate < self.DMSP_ERA:
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
    power_task = HIIPower(**vars(options))
    power_task.run()
