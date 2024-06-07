import argparse
import ee
from datetime import date
from task_base import HIITask
from input_preprocess import HIIPowerPreprocessTask


class HIIPower(HIITask):
    inputs = {
        "nightlights": {
            "ee_type": HIITask.IMAGECOLLECTION,
            "ee_path": "nightlights_path",
            "maxage": 1,
        }
    }

    quantiles = {
        "0": {"value": 0, "min": 0, "max": 9},
        "1": {"value": 1, "min": 10, "max": 10},
        "2": {"value": 2, "min": 11, "max": 12},
        "3": {"value": 3, "min": 13, "max": 13},
        "4": {"value": 4, "min": 14, "max": 16},
        "5": {"value": 5, "min": 17, "max": 19},
        "6": {"value": 6, "min": 20, "max": 25},
        "7": {"value": 7, "min": 26, "max": 35},
        "8": {"value": 8, "min": 36, "max": 51},
        "9": {"value": 9, "min": 52, "max": 62},
        "10": {"value": 10, "min": 63, "max": 63},
    }

    thresholds = {
        "global": 10,
        "latitude": {"min_lat": 45, "min_val": 10, "max_lat": 60, "max_val": 30},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nightlights, _ = self.get_most_recent_image(
            ee.ImageCollection(self.inputs["nightlights"]["ee_path"])
        )
        self.quantiles = ee.Dictionary(self.quantiles)

    def nightlights_path(self):
        if self.taskdate.year < 2013:
            return "projects/sat-io/open-datasets/Harmonized_NTL/dmsp"
        else:
            return "projects/sat-io/open-datasets/Harmonized_NTL/viirs"

    def latitude_mask(self):
        min_lat = self.thresholds["latitude"]["min_lat"] ** 4
        max_lat = self.thresholds["latitude"]["max_lat"] ** 4
        min_val = self.thresholds["latitude"]["min_val"]
        max_val = self.thresholds["latitude"]["max_val"]

        latitude_image = (
            ee.Image.pixelLonLat()
            .select("latitude")
            .abs()
            .pow(4)
            .clamp(min_lat, max_lat)
            .divide(max_lat)
        )
        latitude_threshold_image = (
            latitude_image.subtract(min_lat / max_lat)
            .divide(1 - min_lat / max_lat)
            .multiply(max_val - min_val)
            .add(min_val)
            .rename("latitude_thresholds")
        )
        return latitude_threshold_image

    def calc(self):
        if self.taskdate.year < 2013:
            self.nightlights = self.nightlights.updateMask(
                self.nightlights.gte(self.thresholds["global"])
            )

        else:
            self.nightlights = self.nightlights.updateMask(
                self.nightlights.gte(self.latitude_mask())
            )

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

        self.export_image_ee(hii_power_driver, "driver/power_alt_nl")

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
