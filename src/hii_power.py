import argparse 
import ee
from datetime import datetime, timezone
from task_base import EETask


class HIIPower(EETask):
    ee_rootdir = "projects/HII/v1/sumatra_poc"
    ee_driverdir = 'driver/power'
    # if input lives in ee, it should have an "ee_path" pointing to an ImageCollection/FeatureCollection
    inputs = {
        "viirs": {
            "ee_type": EETask.IMAGECOLLECTION,
            "ee_path": 'NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG',
        },
        "jrc": {
            "ee_type": EETask.IMAGE,
            "ee_path": "JRC/GSW1_0/GlobalSurfaceWater"
        },
        "caspian": {
            "ee_type": EETask.IMAGE,
            "ee_path": "users/aduncan/caspian"
        },
        "ocean": {
            "ee_type": EETask.IMAGE,
            "ee_path": "users/aduncan/cci/ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v40",
        },
        "watermask": {
            "ee_type": EETask.IMAGE,
            "ee_path": "projects/HII/v1/source/phys/watermask_jrc70_cciocean",
        }
            }
    
    gpw_cadence = 5

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_aoi_from_ee("{}/sumatra_poc_aoi".format(self.ee_rootdir))

    def calc(self):

        watermask = ee.Image(self.inputs['watermask']['ee_path'])
        viirs = ee.ImageCollection(self.inputs['viirs']['ee_path'])
        polar = ee.Feature(ee.FeatureCollection('projects/HII/v1/misc/polar_for_viirs_correction').first())

        viirs_polar_correction = ee.FeatureCollection([polar.set({'foo':1})]).reduceToImage(['foo'], ee.Reducer.sum())

        def viirs_func(year):
            startdate = ee.Date.fromYMD(year,1,1)
            viirs_median = viirs.filterDate(startdate,startdate.advance(1,'year'))\
                                          .select('avg_rad')\
                                          .median()\
                                          .reproject(crs='EPSG:4326',scale=300)
                                                         
                                                      
            viirs_outsideSolar = viirs_median.multiply(viirs_polar_correction.neq(1))
            viirs_insideSolar = viirs_median.multiply(viirs_polar_correction.eq(1))
            viirs_noSolar = viirs_outsideSolar.add(viirs_insideSolar)
            return viirs_noSolar

        def calibration(image):
            image = image.unmask();
            return ee.Image(1569.44).add(ee.Image(828.19).multiply(image.log()))

        def reproject_dmsp(image):
            return image.reproject(crs='EPSG:4326',scale=300)

 

        dmsp_1992 = ee.Image('users/aduncan/yale/F101992')
        dmsp_1993 = ee.Image('users/aduncan/yale/F101993')
        dmsp_1994 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F101994'),\
                          ee.Image('users/aduncan/yale/F121994')])\
                            .mean()
        dmsp_1995 = ee.Image('users/aduncan/yale/F121995')
        dmsp_1996 = ee.Image('users/aduncan/yale/F121996')
        dmsp_1997 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F121997'),\
                          ee.Image('users/aduncan/yale/F141997')])\
                            .mean()
        dmsp_1998 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F121998'),\
                          ee.Image('users/aduncan/yale/F141998')])\
                            .mean()
        dmsp_1999 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F121999'),\
                          ee.Image('users/aduncan/yale/F141999')])\
                            .mean()
        dmsp_2000 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F142000'),\
                          ee.Image('users/aduncan/yale/F152000')])\
                            .mean()
        dmsp_2001 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F142001'),\
                          ee.Image('users/aduncan/yale/F152001')])\
                            .mean()
        dmsp_2002 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F142002'),\
                          ee.Image('users/aduncan/yale/F152002')])\
                            .mean()
        dmsp_2003 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F142003'),\
                          ee.Image('users/aduncan/yale/F152003')])\
                            .mean()
        dmsp_2004 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F152004'),\
                          ee.Image('users/aduncan/yale/F162004')])\
                            .mean()
        dmsp_2005 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F152005'),\
                          ee.Image('users/aduncan/yale/F162005')])\
                            .mean()
        dmsp_2006 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F152006'),\
                          ee.Image('users/aduncan/yale/F162006')])\
                            .mean()
        dmsp_2007 = ee.ImageCollection([\
                          ee.Image('users/aduncan/yale/F152007'),\
                          ee.Image('users/aduncan/yale/F162007')])\
                            .mean()
        dmsp_2008 = ee.Image('users/aduncan/yale/F162008')
        dmsp_2009 = ee.Image('users/aduncan/yale/F162009')
        dmsp_2009 = ee.Image('users/aduncan/yale/F162009')
        dmsp_2010 = ee.Image('users/aduncan/yale/F182010')
        dmsp_2011 = ee.Image('users/aduncan/yale/F182011')
        dmsp_2012 = ee.Image('users/aduncan/yale/F182012')






        yearlist = ee.List.sequence(2014,2019)


        '''
        viirs_ic = ee.ImageCollection(yearlist.map(function(year){return viirs_func(year)}))
          .map(function(image){return calibration(image)});
        '''
        viirs_ic = ee.ImageCollection(yearlist.map(viirs_func)).map(calibration)

        # repeating for 2013... could interpolate later?
        dmsp_ic = ee.ImageCollection([\
                  dmsp_1992,\
                  dmsp_1993,\
                  dmsp_1994,\
                  dmsp_1995,\
                  dmsp_1996,\
                  dmsp_1997,\
                  dmsp_1998,\
                  dmsp_1999,\
                  dmsp_2000,\
                  dmsp_2001,\
                  dmsp_2002,\
                  dmsp_2003,\
                  dmsp_2004,\
                  dmsp_2005,\
                  dmsp_2006,\
                  dmsp_2007,\
                  dmsp_2008,\
                  dmsp_2009,\
                  dmsp_2010,\
                  dmsp_2011,\
                  dmsp_2012,\
                  dmsp_2012\
                    ]).map(reproject_dmsp)
  
        # NEED TO SORT DESCENDING!!! This was a map but don't need right now...
        ee_taskdate = ee.Date(self.taskdate.strftime(self.DATE_FORMAT))
        year_num = ee.Number(ee_taskdate.get('year')).subtract(1992)

        hii_power_oneyear = ee.Image(dmsp_ic.merge(viirs_ic).toList(30).get(year_num))

        hii_power_driver = hii_power_oneyear.unmask(0).multiply(0.01).updateMask(watermask).multiply(10)
        self.export_image_ee(hii_power_driver, '{}/{}'.format(self.ee_driverdir, 'hii_power_driver'))

    def check_inputs(self):
        super().check_inputs()
        # add any task-specific checks here, and set self.status = self.FAILED if any fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--taskdate', default=datetime.now(timezone.utc).date())
    options = parser.parse_args()
    power_task = HIIPower(**vars(options))
    power_task.run()
