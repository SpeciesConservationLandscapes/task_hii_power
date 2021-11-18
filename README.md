HII POWER DRIVER
---------------

## What does this task do?

This task calculates the (unitless) "influence" of power on the terrestrial surface as one of the key drivers for a combined [Human Influence Index](https://github.com/SpeciesConservationLandscapes/task_hii_weightedsum). "Influence" is a pressure score based on the intensity of electricity usage as measured by "Night Time Lights" datasets. An inter-calibrated nighttime lights dataset with values ranging from 0 - 63 was transformed to a pressure score of 0 - 10 using equal interval quantiles, calculated form the earliest available image (1992). These were:

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

### Input Dataset Callibration
Two distinct nightlight datasets are integrated to produce an annually calibrated dataset of nightlights from 1992 through the present. These the [Defense Meteorological Satellite Program (DMSP)/Operational Linescan System (OLS)](https://eogdata.mines.edu/products/dmsp/) and the [Visible Infrared Imaging Radiometry Suite (VIIRS) on the Suomi National Polar-orbiting Partnership Satellite](https://eogdata.mines.edu/products/vnl/). DMSP provides data from 1992 - 2013 and VIIRS provide data from 2012 through the present.

Inconsistencies within the DMSP time series requires implementing intra-calibration within the DMPS dataset. Here the intra-calibrated dataset produced by [Li et al. 2020](https://www.nature.com/articles/s41597-020-0510-y) is used. Key differences between the DMSP and VIIRS datasets also require inter-calibration of the two datasets. In addition to the calibration process to match VIIRS to DMSP, VIIRS requires substantial noise reduction. These key differences are:

| Description | DMSP | VIIRS | Calibrated |
| :--- | :--- | :--- | :--- |
| Date Range | 1992 - 2013 | 2012 - present | 1992 - present |
| Temporal Resolution | Annual | Monthly | Annual |
| Spatial Resolution | 30 arc-second | 15 arc-seconds | 30 arc-second |
| Pixel Value | Digital Number | Radiance | Digital Number |
| Data Range | 0 - 63 | -1.5 - 193564.92 | 0 - 63 |

 The steps to noise removal and calibration steps are as follows:

 1. Annual composites are created using the median value of all pixels for a given year with raw values greater than 0.

 2. General noise (which generally increases with distance from the equator) is removed by creating an image mask using a pixel latitude image. An exponential function is used to stretch the threshold values from a minimum to maximum threshold of 0.1 and 0.75 respectively from 0° to ±60°. The formula for this is:

&emsp; &emsp; &emsp; t = (|lat|<sup>4</sup> / 60<sup>4</sup>) x (t<sub>max</sub> - t<sub>min</sub>) + (t<sub>min</sub>)

&emsp; &emsp; &emsp; &emsp; Where: <br />
&emsp; &emsp; &emsp; &emsp; t is the final threshold image <br />
&emsp; &emsp; &emsp; &emsp; lat is a latitude image <br />
&emsp; &emsp; &emsp; &emsp; t<sub>min</sub> is the minimum threshold <br />
&emsp; &emsp; &emsp; &emsp; t<sub>max</sub> is the maximum threshold

3. Regression coeffecients were calculated to convert VIIRS values to DMPS values by comparing pixel values of areas with stable lights through time. These were determined by selecting pixels from the DMSP time series from 2000 - 2012 that had a standard deviation less than 2. A stratified sample of 200 pixel values for each available DMSP value was calculated. The resulting equation to transform VIIRS values to match DMSP was determined to be:

&emsp; &emsp; &emsp; calibrated_viirs = log(viirs) x 10.53 + 24.62

4. Values in the calibrated image greater than 63 were set 63.

5. For the year 2012, in which both DMSP and VIIRS images are available, the final calibrated image was calculated as the mean of DMPS and VIIRS for 2012.
