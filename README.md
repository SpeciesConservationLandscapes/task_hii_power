HII POWER DRIVER
---------------

## What does this task do?

This task calculates the impact of power on the terrestrial surface as one of the key drivers for a combined 
[Human Impact Index](https://github.com/SpeciesConservationLandscapes/task_hii_weightedsum). 
"Impact" is a pressure score based on the intensity of electricity usage as measured by "Night Time Lights" datasets. 
An inter-calibrated nighttime lights dataset with values ranging from 0 - 63 is transformed to a pressure score of 
0 - 10 using equal interval quantiles, calculated from the earliest available image (1992) using the `quantile_calc` 
function in `input_preprocess.py`. These are:

 ```
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
 ```
The output HII driver calculated by this task is, like all other HII drivers, unitless; it refers to an absolute 0-10
scale but is not normalized to it, so the actual range of values may be smaller than 0-10.

## Input Dataset Calibration
The HII power driver (0-10) is calculated from the previous calendar year's calibrated nightlights dataset, which is 
produced on demand by the task by using the `HIIPowerPreprocessTask` task. Two distinct source nightlight datasets are 
used to calculate this calibrated version following the methods below; these are the 
Defense Meteorological Satellite Program (DMSP)/Operational Linescan System (OLS) and the 
Visible Infrared Imaging Radiometry Suite (VIIRS) on the Suomi National Polar-orbiting Partnership Satellite. 
DMSP provides data from 1992 - 2013 and VIIRS provides data from 2012 through the present.

Inconsistencies within the DMSP time series require implementing intra-calibration within the DMSP dataset. 
Here the intra-calibrated dataset produced by [Li et al. 2020](https://www.nature.com/articles/s41597-020-0510-y) 
is used. Key differences between the DMSP and VIIRS datasets also require inter-calibration of the two datasets. 
In addition to the calibration process to match VIIRS to DMSP, VIIRS requires substantial noise reduction. 
These key differences are:

| Description | DMSP | VIIRS | Calibrated |
| :--- | :--- | :--- | :--- |
| Date Range | 1992 - 2013 | 2012 - present | 1992 - present |
| Temporal Resolution | Annual | Monthly | Annual |
| Spatial Resolution | 30 arc-second | 15 arc-seconds | 30 arc-second |
| Pixel Value | Digital Number | Radiance | Digital Number |
| Data Range | 0 - 63 | -1.5 - 193564.92 | 0 - 63 |

### Noise removal

1. Annual composites are created using the median value of all pixels for a given year with raw values greater than 0.
2. General noise (which generally increases with distance from the equator) is removed by creating an image mask using 
   a pixel latitude image. An exponential function is used to stretch the threshold values from a minimum to 
   maximum threshold of 0.1 and 0.75 respectively from 0° to ±60°. The formula for this is:

&emsp; &emsp; &emsp; t = (|lat|<sup>4</sup> / 60<sup>4</sup>) x (t<sub>max</sub> - t<sub>min</sub>) + (t<sub>min</sub>)

&emsp; &emsp; &emsp; &emsp; Where: <br />
&emsp; &emsp; &emsp; &emsp; t is the final threshold image <br />
&emsp; &emsp; &emsp; &emsp; lat is a latitude image <br />
&emsp; &emsp; &emsp; &emsp; t<sub>min</sub> is the minimum threshold <br />
&emsp; &emsp; &emsp; &emsp; t<sub>max</sub> is the maximum threshold

### Calibration steps

1. Regression coefficients are calculated to convert VIIRS values (beginning with 2014) to DMSP values by comparing 
   pixel values of areas with stable lights through time. These are determined statically by selecting pixels from 
   the DMSP time series from 2000 - 2013 that have a standard deviation less than 2. A stratified sample of 300 pixel 
   values for each available DMSP value as of the latest available image (2013) is calculated and used to derive 
   constant slope and intercept coefficients. See the `stable_light_points` function in `input_preprocess.py`.
2. The resulting equation to transform VIIRS values to match DMSP is determined to be:  
&emsp; &emsp; &emsp; calibrated_viirs = log(viirs) * 11.62 + 19.4 [R<sup>2</sup>: 67.16]
3. Values in the calibrated image greater than 63 are set 63, matching the range of values from DMSP.

The result and goal is an annual night lights dataset derived from VIIRS that is calibrated to the pre-2013 
DMSP values, which is then transformed into the HII power driver by the main task.

## Variables and Defaults

### Environment variables
```
SERVICE_ACCOUNT_KEY=<GOOGLE SERVICE ACCOUNT KEY>
```

### Class constants

```
scale=300
```

## Usage

*All parameters may be specified in the environment as well as the command line.*

```
/app # python task.py --help
usage: task.py [-h] [-d TASKDATE] [--overwrite]

optional arguments:
  -h, --help            show this help message and exit
  -d TASKDATE, --taskdate TASKDATE
  --overwrite           overwrite existing outputs instead of incrementing
```

In addition, the `HIIPowerPreprocessTask` task can be called to run one of three preprocessing steps used in the 
main `task.py`, with the following constants available as arguments:
```
BUCKET = "hii-scratch"
CALC_PREVIOUS_ANNUAL_VIIRS = "viirs"
CALC_CALIBRATION_COEFFICIENTS = "coefficients"
CALC_WEIGHTING_QUANTILES = "quantiles"
```
The first of these, `viirs`, is used on demand by the main power `task.py`.
```
/app# python input_preprocess.py --help
usage: input_preprocess.py [-h] [-d TASKDATE] [-j [{viirs,coefficients,quantiles}]] [--overwrite]

optional arguments:
  -h, --help            show this help message and exit
  -d TASKDATE, --taskdate TASKDATE
  -j [{viirs,coefficients,quantiles}], --job [{viirs,coefficients,quantiles}]
                        Calculation to run with HIIPowerPreprocessTask.
  --overwrite           overwrite existing outputs instead of incrementing
```

### License
Copyright (C) 2022 Wildlife Conservation Society
The files in this repository  are part of the task framework for calculating 
Human Impact Index and Species Conservation Landscapes (https://github.com/SpeciesConservationLandscapes) 
and are released under the GPL license:
https://www.gnu.org/licenses/#GPL
See [LICENSE](./LICENSE) for details.
