import argparse
import ee
from task_base import HIITask


class HIIPower(HIITask):
    scale = 300

    inputs = {
        "dmsp_viirs_calibrated": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "projects/HII/v1/source/nightlights/dmsp_viirs_calibrated",
            "maxage": 1,
        },
        "watermask": {
            "ee_type": HIITask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
            "static": True,
        },
    }
    quantiles = {
        "0": {"value": 0, "min": 0, "max": 0},
        "1": {"value": 1, "min": 1, "max": 4},
        "2": {"value": 2, "min": 5, "max": 5},
        "3": {"value": 3, "min": 6, "max": 6},
        "4": {"value": 4, "min": 7, "max": 7},
        "5": {"value": 5, "min": 8, "max": 9},
        "6": {"value": 6, "min": 10, "max": 11},
        "7": {"value": 7, "min": 12, "max": 16},
        "8": {"value": 8, "min": 17, "max": 30},
        "9": {"value": 9, "min": 31, "max": 62},
        "10": {"value": 10, "min": 63, "max": 63},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nightlights, _ = self.get_most_recent_image(
            ee.ImageCollection(self.inputs["dmsp_viirs_calibrated"]["ee_path"])
        )
        self.watermask = ee.Image(self.inputs["watermask"]["ee_path"])
        self.quantiles = ee.Dictionary(self.quantiles)

    # TODO: add code for calculating regression coeffecients used in callibration
    # TODO: add code for exporting callibrated nightlights, including noise correction for 2013 - current
    # TODO: include method for checking if nightlight image is available and export if not

    def calc(self):
        def power_classify(value):
            bin = ee.Dictionary(self.quantiles.get(value))
            power_value = ee.Number(bin.get("value"))
            min = ee.Number(bin.get("min"))
            max = ee.Number(bin.get("max"))
            return (
                self.nightlights.gte(min)
                .And(self.nightlights.lte(max))
                .selfMask()
                .multiply(power_value)
                .unmask(0)
                .int()
            )

        hii_power_driver = (
            ee.ImageCollection(self.quantiles.keys().map(power_classify))
            .reduce(ee.Reducer.max())
            .updateMask(self.watermask)
            .multiply(100)
            .int()
            .rename("hii_power_driver")
        )

        self.export_image_ee(hii_power_driver, "driver/power")

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
    power_task = HIIPower(**vars(options))
    power_task.run()
