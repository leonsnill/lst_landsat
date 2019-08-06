# Landsat Land Surface Temperature (LST) Retrieval – Google Earth Engine (GEE) Implementation

Method used by ...(2019) (link/doi)

Google Earth Engine (GEE) Implementation of the Single-Channel (SC) algorithm developed by Jiménez-Muñoz & Sobrino (2003), 
Jiménez-Muñoz et al. (2009) and Jiménez-Muñoz et al. (2014) for retrieving statistical metrics of LST
(and multispectral indices) from Landsat TM, ETM+ and OLI-TIRS data.
The atmospheric functions used in the LST algorithm are approximated using data on atmospheric water vapor content from the NCEP/NCAR Reanalysis Data available within the GEE data catalogue. Currently, the approximation is optimized for high-latitude regions with usually low water vapor content.
Surface emissivity is calculated using the Simplified Normalized Difference Vegetation Index Threshold (SNDVI) approach as described by Sobrino et al. (2008).

### Prerequisites

A succsessful installation of the Earth Engine Python API is required (see https://developers.google.com/earth-engine/python_install for instructions).
After the installation, the following commands should be executable without any error messages:

```
import ee
ee.Initialize()
```

### Functionality

The script currently allows for retrieving temporal statistical metrics (e.g. mean, stdDev, ...) within defined space and time. It encompasses pre-processing the Landsat data including cloud, cloud shadow and snow masking using the internal quality assessment band from the CFMask algorithm of the surface reflectance products.
The user defined surface parameters (list down below) will then be retrieved and the statistical metrics calculated for each of them.

#### Surface Parameters and Statistical Metrics

Currently the following parameters and metrics are implemented (with the possibility to be easily extended):
- Land Surface Temperature ('LST'), Normalized Difference Vegetation Index ('NDVI'), Normalized Difference Water Index after Gao (1996) ('NDWI'), Tasseled Cap Greenness ('TCG'), Brightness ('TCB') and Wetness ('TCW') using the surface reflectance coefficients from Christ (1985).

Example, how to add another index to the script (here: Enhanced Vegetation Index; EVI):

Define function in 'FUNCTIONS' section:
```
def fun_evi(img):
    evi = img.expression(
                         'Gain*((NIR-R)/(NIR + C1*R - C2*B + L))',
                         {
                         'Gain': ee.Number(2.5),
                         'B': img.select(['B']),
                         'R': img.select(['R']),
                         'NIR': img.select(['NIR']),
                         'L': ee.Number(1),
                         'C1': ee.Number(6),
                         'C2': ee.Number(7.5)
                         }).rename('EVI')
    return img.addBands(evi)
```
Add calculation to other index calculations:
```
imgCol_merge = imgCol_merge.map(fun_evi)
```
Lastly, be aware of the value range of the index (e.g. -1 to 1 for the EVI) and if needed edit the following section to ensure the output is properly scaled and defined:
```
    for metric in select_metrics:
        if metric == 'ts':
            temp = imgCol_mosaic.select(['TIME', parameter]).reduce(ee.Reducer.sensSlope())
            temp = temp.select('slope')
            temp = temp.multiply(365)
            temp = temp.multiply(100000000).int32()
        elif metric == 'nobs':
            temp = imgCol_mosaic.select(parameter).count()
            temp = temp.int16()
        else:
            if metric == 'percentile':
                temp = imgCol_mosaic.select(parameter).reduce(ee.Reducer.percentile(percentiles))
            else:
                reducer = lookup_metrics[metric]
                temp = imgCol_mosaic.select(parameter).reduce(reducer)
            if parameter == 'LST':
                temp = temp.multiply(100).int16()
            else:
                temp = temp.multiply(10000).int16()


```
In this example, nothing needs to be changed in the above. Finally, one can retrieve the newly defined parameter by simply including the band name (here: 'EVI') in the 'select_parameters'-list, e.g.:
```
select_parameters = ['TCG', 'TCB', 'TCW', 'EVI']
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
