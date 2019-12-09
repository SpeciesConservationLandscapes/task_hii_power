import ee
from task_base import EETask


class HIILanduse(EETask):
    ee_rootdir = "projects/HII/v1/sumatra_poc"
    ee_driverdir = 'driver/landuse'
    # if input lives in ee, it should have an "ee_path" pointing to an ImageCollection/FeatureCollection
    inputs = {
        "gpw": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": "CIESIN/GPWv411/GPW_Population_Density",
        },
        "jrc": {
            "ee_type": EETask.IMAGE,
            "ee_path": "JRC/GSW1_0/GlobalSurfaceWater"
        },
        "caspian": {
            "ee_type": EETask.IMAGE,
            "ee_path": "users/aduncan/caspian"
        },
        "gpw_2015": {
            "ee_type": EETask.IMAGE,
            "ee_path": "CIESIN/GPWv411/GPW_Population_Density/gpw_v4_population_density_rev11_2015_30_sec"
        },
        "esacci": {
            "ee_type": EETask.IMAGE,
            "ee_path": "users/aduncan/cci/ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v207"
        }
    }
    

    def FUNC_IC(item,image_coll):
        year_property = ee.Number(item).add(1992)
        single_year = ee.ImageCollection([ESACCI.select([item]).set('year',year_property)])
        return ee.ImageCollection(image_coll).merge(single_year)
  



    gpw_cadence = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))

    #def test_function(number):
    #    print(number)

    #print(test_function(5555555555555))

    def square(x):
      square = x * x
      return square

    test1 = 'test1'

    def something(self):
      something_else = self.test1
      return something_else

    def calc(self):
        print(self.test1) # works
        print(self.something()) # works

        ESACCI = ee.ImageCollection(self.inputs['esacci']['ee_path'])
        year_list =   ee.List.sequence(0,23)
        ESACCI_ic =   ee.ImageCollection(year_list.iterate(FUNC_IC(), ee.ImageCollection([])))

        print(year_list.getInfo())





        gpw_2015 = ee.ImageCollection(self.inputs['gpw_2015']['ee_path'])

        ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        gpw_prior = gpw.filterDate(ee_taskdate.advance(-self.gpw_cadence, 'year'), ee_taskdate).first()
        gpw_later = gpw.filterDate(ee_taskdate, ee_taskdate.advance(self.gpw_cadence, 'year')).first()
        gpw_diff = gpw_later.subtract(gpw_prior)
        numerator = ee_taskdate.difference(gpw_prior.date(), 'day')
        gpw_diff_fraction = gpw_diff.multiply(numerator.divide(self.gpw_cadence * 365))
        gpw_taskdate = gpw_prior.add(gpw_diff_fraction)
        gpw_taskdate_300m = gpw_taskdate.resample().reproject(crs=self.crs, scale=self.scale)

        gpw_venter = gpw_taskdate_300m.add(ee.Image(1))\
            .log()\
            .multiply(ee.Image(3.333))
        # TODO: mask water with centralized HII-defined water images

        #self.export_image_ee(gpw_venter, '{}/{}'.format(self.ee_driverdir, 'hii_popdens_driver'))

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    landuse_task = HIILanduse()
    landuse_task.run()
