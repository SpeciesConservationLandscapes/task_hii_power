import argparse
import ee
from datetime import datetime, timezone
from task_base import HIITask


class HIIPower(HIITask):
    ee_rootdir = "projects/HII/v1"
    ee_driverdir = "driver/power"
    inputs = {
        "watermask": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True,
        },
        "dmsp_viirs_merged": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/source/nightlights/dmsp_viirs_merged",
            "maxage": 3,
        },
    }
    scale = 300
    gpw_cadence = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.realm = kwargs.pop("realm", None)
        self.set_aoi_from_ee("projects/HII/v1/source/realms/" + self.realm)

    def calc(self):
        nightlights, nightlights_date = self.get_most_recent_image(
            ee.ImageCollection(self.inputs["dmsp_viirs_merged"]["ee_path"])
        )
        watermask = ee.Image(self.inputs["watermask"]["ee_path"])

        # Adam TODO question: Should I just find the julian day, divide by 365, and use that to find a date's value between two 01-01s?
        hii_power_driver = nightlights.multiply(0.01).updateMask(watermask)
        self.export_image_ee(
            hii_power_driver, "{}/{}".format(self.ee_driverdir, "aois/" + self.realm)
        )

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--realm", default="Afrotropic")
    parser.add_argument("-d", "--taskdate")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing outputs instead of incrementing",
    )
    options = parser.parse_args()
    power_task = HIIPower(**vars(options))
    power_task.run()
