'''
LANDSAT LAND SURFACE TEMPERATURE AND MULTISPECTRAL INDICES METRICS RETRIEVAL
Leon Nill (2019)

GEE Implementation of the Single-Channel (SC) algorithm developed by Jiménez-Muñoz & Sobrino (2003),
Jiménez-Muñoz et al. (2009) and Jiménez-Muñoz et al. (2014) for retrieving statistical metrics of LST
from Landsat TM, ETM+ and OLI-TIRS data. The atmospheric functions used in the algorithm are
approximated using data on atmospheric water vapor content from the NCEP/NCAR Reanalysis Data.
Currently, the approximation is optimized for high-latitude regions with usually low water vapor
content. Surface emissivity is calculated using the Simplified Normalized Difference Vegetation Index
Threshold (SNDVI) approach as described by Sobrino et al. (2008).

  User Requirements:
  select_parameter          [STRING] List of selected parameter to process
                            Currently implemented: 'LST', 'NDVI', 'NDWI', 'TCG', 'TCB', 'TCW'
                            Output scaled by factor 10,000, except LST scaled by 0.2 and 1000 both in Int16 GTiff.
  select_metrics            [STRING] List of selected metrics to calculate
                            Currently implemented: 'mean', 'median', 'min', 'max', 'std', 'percentile',
                            'ts' (Theil-Sen slope), 'nobs' (Number of Observations).
  percentiles               [INT] List of percentiles to calculate if 'percentile' specified in 'select_metrics'.
  year_start, year_end      [INT] Start- and end-year of the temporal integration.
  month_start, month_end    [INT] Start- and end-month of the temporal integration.
  select_roi                [xMin, yMin, xMax, yMax] List of corner coordinates.
  max_cloud_cover           [INT/FLOAT] Maximum scene (land) cloud cover to be included in selection.
  epsg                      [STRING] EPSG code for the target spatial reference, e.g. 'EPSG:4326'.
  pixel_resolution          [INT/FLOAT] Output pixel resolution in meters.
  roi_filename              [STRING] Name to be included in output filename, e.g. abbreviation of study area.
  filename                  [STRING] Output master-filename. Specific metrics are automatically appended.


  Algorithm Specifications:
  ndvi_v, ndvi_s            [FLOAT] Maximum/minimum NDVI values at which full vegetation (_v)/soil (_s)
                            cover exists. I.e. maximum/minumum emissivity is reached.
  epsilon_v, epsilon_s,     [FLOAT] Emissivity values of full vegetation (_v), soil (_s) or water
  epsilon_w                 (_w) cover.
  t_threshold               [INT/FLOAT] Minimum temperature in degrees Celcius below which values are masked.
                            Compensates for issues of insufficient cloud detection in CFMASK.
  cs_*                      [FLOAT] List of coefficients used to approximate water vapor content.

'''


import ee
ee.Initialize()


# ====================================================================================================#
# INPUT
# ====================================================================================================#

# --------------------------------------------------
# User Requirements
# --------------------------------------------------
select_parameters = ['LST', 'NDWI']
select_metrics = ['mean', 'nobs', 'ts']
percentiles = [10, 90]

# Time
year_start = 1985
year_end = 2018
month_start = 7
month_end = 8

# Space
select_roi = ee.Geometry.Rectangle([-137.67, 67.40, -132.21, 69.84])
max_cloud_cover = 60
epsg = 'EPSG:32608'
pixel_resolution = 30
roi_filename = 'MRD'

# --------------------------------------------------
# Algorithm Specifications
# --------------------------------------------------
ndvi_v = 0.6
ndvi_s = 0.2

epsilon_v = 0.985
epsilon_s = 0.97
epsilon_w = 0.99

t_threshold = 8

# Jiménez‐Muñoz et al. (2009) (TM & ETM+) TIGR1761 and Jiménez‐Muñoz et al. (2014) OLI-TIRS GAPRI4838
cs_l8 = [0.04019, 0.02916, 1.01523,
         -0.38333, -1.50294, 0.20324,
         0.00918, 1.36072, -0.27514]
cs_l7 = [0.06518, 0.00683, 1.02717,
         -0.53003, -1.25866, 0.10490,
         -0.01965, 1.36947, -0.24310]
cs_l5 = [0.07518, -0.00492, 1.03189,
         -0.59600, -1.22554, 0.08104,
         -0.02767, 1.43740, -0.25844]

# ====================================================================================================#
# FUNCTIONS
# ====================================================================================================#

lookup_metrics = {
    'mean': ee.Reducer.mean(),
    'min': ee.Reducer.min(),
    'max': ee.Reducer.max(),
    'std': ee.Reducer.stdDev(),
    'median': ee.Reducer.median(),
    'ts': ee.Reducer.sensSlope()
}


# --------------------------------------------------
# RENAME BANDS
# --------------------------------------------------

def fun_bands_l57(img):
       bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
       thermal_band = ['B6']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       new_thermal_bands = ['TIR']
       vnirswir = img.select(bands).multiply(0.0001).rename(new_bands)
       tir = img.select(thermal_band).multiply(0.1).rename(new_thermal_bands)
       return vnirswir.addBands(tir).copyProperties(img, ['system:time_start'])


def fun_bands_l8(img):
       bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
       thermal_band = ['B10']
       new_bands = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']
       new_thermal_bands = ['TIR']
       vnirswir = img.select(bands).multiply(0.0001).rename(new_bands)
       tir = img.select(thermal_band).multiply(0.1).rename(new_thermal_bands)
       return vnirswir.addBands(tir).copyProperties(img, ['system:time_start'])


# --------------------------------------------------
# MASKING
# --------------------------------------------------

# Function to cloud mask Landsat TM, ETM+, OLI_TIRS Surface Reflectance Products
def fun_mask_ls_sr(img):
       cloudShadowBitMask = ee.Number(2).pow(3).int()
       cloudsBitMask = ee.Number(2).pow(5).int()
       snowBitMask = ee.Number(2).pow(4).int()
       qa = img.select('pixel_qa')
       mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(
              qa.bitwiseAnd(cloudsBitMask).eq(0)).And(
              qa.bitwiseAnd(snowBitMask).eq(0))
       return img.updateMask(mask)


# Function to mask LST below certain temperature threshold
def fun_mask_T(img):
    mask = img.select('LST').gt(t_threshold)
    return img.updateMask(mask)



# --------------------------------------------------
# MATCHING AND CALIBRATION
# --------------------------------------------------

# Radiometric Calibration
def fun_radcal(img):
    radiance = ee.Algorithms.Landsat.calibratedRadiance(img).rename('RADIANCE')
    return img.addBands(radiance)


# L to ee.Image
def fun_l_addband(img):
    l = ee.Image(img.get('L')).select('RADIANCE').rename('L')
    return img.addBands(l)


# Create maxDifference-filter to match TOA and SR products
maxDiffFilter = ee.Filter.maxDifference(
    difference=2 * 24 * 60 * 60 * 1000,
    leftField= 'system:time_start',
    rightField= 'system:time_start'
)

# Define join: Water vapor
join_wv = ee.Join.saveBest(
    matchKey = 'WV',
    measureKey = 'timeDiff'
)

# Define join: Radiance
join_l = ee.Join.saveBest(
    matchKey = 'L',
    measureKey = 'timeDiff'
)


# --------------------------------------------------
# PARAMETER CALCULATION
# --------------------------------------------------

# NDVI
def fun_ndvi(img):
    ndvi = img.normalizedDifference(['NIR', 'R']).rename('NDVI')
    return img.addBands(ndvi)


def fun_ndwi(img):
    ndwi = img.normalizedDifference(['NIR', 'SWIR1']).rename('NDWI')
    return img.addBands(ndwi)


# Tasseled Cap Transformation (brightness, greenness, wetness) based on Christ (1985)
def fun_tcg(img):
    tcg = img.expression(
                         'B*(-0.1603) + G*(-0.2819) + R*(-0.4934) + NIR*0.7940 + SWIR1*(-0.0002) + SWIR2*(-0.1446)',
                         {
                         'B': img.select(['B']),
                         'G': img.select(['G']),
                         'R': img.select(['R']),
                         'NIR': img.select(['NIR']),
                         'SWIR1': img.select(['SWIR1']),
                         'SWIR2': img.select(['SWIR2'])
                         }).rename('TCG')
    return img.addBands(tcg)


def fun_tcb(img):
    tcb = img.expression(
                         'B*0.2043 + G*0.4158 + R*0.5524 + NIR*0.5741 + SWIR1*0.3124 + SWIR2*0.2303',
                         {
                         'B': img.select(['B']),
                         'G': img.select(['G']),
                         'R': img.select(['R']),
                         'NIR': img.select(['NIR']),
                         'SWIR1': img.select(['SWIR1']),
                         'SWIR2': img.select(['SWIR2'])
                         }).rename('TCB')
    return img.addBands(tcb)


def fun_tcw(img):
       tcw = img.expression(
              'B*0.0315 + G*0.2021 + R*0.3102 + NIR*0.1594 + SWIR1*(-0.6806) + SWIR2*(-0.6109)',
              {
                     'B': img.select(['B']),
                     'G': img.select(['G']),
                     'R': img.select(['R']),
                     'NIR': img.select(['NIR']),
                     'SWIR1': img.select(['SWIR1']),
                     'SWIR2': img.select(['SWIR2'])
              }).rename('TCW')
       return img.addBands(tcw)


# Fraction Vegetation Cover (FVC)
def fun_fvc(img):
    fvc = img.expression(
        '((NDVI-NDVI_s)/(NDVI_v-NDVI_s))**2',
        {
            'NDVI': img.select('NDVI'),
            'NDVI_s': ndvi_s,
            'NDVI_v': ndvi_v
        }
    ).rename('FVC')
    return img.addBands(fvc)


# Scale Emissivity (Epsilon) between NDVI_s and NDVI_v
def fun_epsilon_scale(img):
    epsilon_scale = img.expression(
        'epsilon_s+(epsilon_v-epsilon_s)*FVC',
        {
            'FVC': img.select('FVC'),
            'epsilon_s': epsilon_s,
            'epsilon_v': epsilon_v
        }
    ).rename('EPSILON_SCALE')
    return img.addBands(epsilon_scale)


# Emissivity (Epsilon)
def fun_epsilon(img):
    pseudo = img.select(['NDVI']).set('system:time_start', img.get('system:time_start'))
    epsilon = pseudo.where(img.expression('NDVI > NDVI_v',
                                          {'NDVI': img.select('NDVI'),
                                           'NDVI_v': ndvi_v}), epsilon_v)
    epsilon = epsilon.where(img.expression('NDVI < NDVI_s && NDVI >= 0',
                                           {'NDVI': img.select('NDVI'),
                                            'NDVI_s': ndvi_s}), epsilon_s)
    epsilon = epsilon.where(img.expression('NDVI < 0',
                                           {'NDVI': img.select('NDVI')}), epsilon_w)
    epsilon = epsilon.where(img.expression('NDVI <= NDVI_v && NDVI >= NDVI_s',
                                           {'NDVI': img.select('NDVI'),
                                            'NDVI_v': ndvi_v,
                                            'NDVI_s': ndvi_s}), img.select('EPSILON_SCALE')).rename('EPSILON')
    return img.addBands(epsilon)


# Function to scale WV content product
def fun_wv_scale(img):
    wv_scaled = ee.Image(img.get('WV')).multiply(0.1).rename('WV_SCALED')
    wv_scaled = wv_scaled.resample('bilinear')
    return img.addBands(wv_scaled)


# --------------------------------------------------
# LAND SURFACE TEMPERATURE CALCULATION
# --------------------------------------------------

# Atmospheric Functions
def fun_af1(cs):
    def wrap(img):
        af1 = img.expression(
            '('+str(cs[0])+'*(WV**2))+('+str(cs[1])+'*WV)+('+str(cs[2])+')',
            {
                'WV': img.select('WV_SCALED')
            }
        ).rename('AF1')
        return img.addBands(af1)
    return wrap


def fun_af2(cs):
    def wrap(img):
        af2 = img.expression(
            '('+str(cs[3])+'*(WV**2))+('+str(cs[4])+'*WV)+('+str(cs[5])+')',
            {
                'WV': img.select('WV_SCALED')
            }
        ).rename('AF2')
        return img.addBands(af2)
    return wrap


def fun_af3(cs):
    def wrap(img):
        af3 = img.expression(
            '('+str(cs[6])+'*(WV**2))+('+str(cs[7])+'*WV)+('+str(cs[8])+')',
            {
                'WV': img.select('WV_SCALED')
            }
        ).rename('AF3')
        return img.addBands(af3)
    return wrap


# Gamma Functions
def fun_gamma_l8(img):
    gamma = img.expression('(BT**2)/(1324*L)',
                           {'BT': img.select('TIR'),
                            'L': img.select('L')
                            }).rename('GAMMA')
    return img.addBands(gamma)


def fun_gamma_l7(img):
    gamma = img.expression('(BT**2)/(1277*L)',
                           {'BT': img.select('TIR'),
                            'L': img.select('L')
                            }).rename('GAMMA')
    return img.addBands(gamma)


def fun_gamma_l5(img):
    gamma = img.expression('(BT**2)/(1256*L)',
                           {'BT': img.select('TIR'),
                            'L': img.select('L')
                            }).rename('GAMMA')
    return img.addBands(gamma)


# Delta Functions
def fun_delta_l8(img):
    delta = img.expression('BT-((BT**2)/1324)',
                           {'BT': img.select('TIR')
                            }).rename('DELTA')
    return img.addBands(delta)


def fun_delta_l7(img):
    delta = img.expression('BT-((BT**2)/1277)',
                           {'BT': img.select('TIR')
                            }).rename('DELTA')
    return img.addBands(delta)


def fun_delta_l5(img):
    delta = img.expression('BT-((BT**2)/1256)',
                           {'BT': img.select('TIR')
                            }).rename('DELTA')
    return img.addBands(delta)


# Land Surface Temperature
def fun_lst(img):
    lst = img.expression(
        '(GAMMA*(((1/EPSILON)*(AF1*L+AF2))+AF3)+DELTA)-273.15',
        {
            'GAMMA': img.select('GAMMA'),
            'DELTA': img.select('DELTA'),
            'EPSILON': img.select('EPSILON'),
            'AF1': img.select('AF1'),
            'AF2': img.select('AF2'),
            'AF3': img.select('AF3'),
            'L': img.select('L')
        }
    ).rename('LST')
    return img.addBands(lst)


def fun_mask_lst(img):
    mask = img.select('LST').gt(t_threshold)
    return img.updateMask(mask)


# --------------------------------------------------
# MOSAICKING
# --------------------------------------------------
def fun_date(img):
    return ee.Date(ee.Image(img).date().format("YYYY-MM-dd"))


def fun_getdates(imgCol):
    return ee.List(imgCol.toList(imgCol.size()).map(fun_date))


def fun_mosaic(date, newList):
    # cast list & date
    newList = ee.List(newList)
    date = ee.Date(date)

    # filter img-collection
    filtered = ee.ImageCollection(subCol.filterDate(date, date.advance(1, 'day')))

    # check duplicate
    img_previous = ee.Image(newList.get(-1))
    img_previous_datestring = img_previous.date().format("YYYY-MM-dd")
    img_previous_millis = ee.Number(ee.Date(img_previous_datestring).millis())

    img_new_datestring = filtered.select(parameter).first().date().format("YYYY-MM-dd")
    img_new_date = ee.Date(img_new_datestring).millis()
    img_new_millis = ee.Number(ee.Date(img_new_datestring).millis())

    date_diff = img_previous_millis.subtract(img_new_millis)

    # mosaic
    img_mosaic = ee.Algorithms.If(
        date_diff.neq(0),
        filtered.select(parameter).mosaic().set('system:time_start', img_new_date),
        ee.ImageCollection(subCol.filterDate(pseudodate, pseudodate.advance(1, 'day')))
    )

    tester = ee.Algorithms.If(date_diff.neq(0), ee.Number(1), ee.Number(0))

    return ee.Algorithms.If(tester, newList.add(img_mosaic), newList)


def fun_timeband(img):
    time = ee.Image(img.metadata('system:time_start', 'TIME').divide(86400000))
    timeband = time.updateMask(img.select(parameter).mask())
    return img.addBands(timeband)


# ====================================================================================================#
# EXECUTE
# ====================================================================================================#


# --------------------------------------------------
# TRANSFORM CLIENT TO SERVER SIDE
# --------------------------------------------------
ndvi_v = ee.Number(ndvi_v)
ndvi_s = ee.Number(ndvi_s)

epsilon_v = ee.Number(epsilon_v)
epsilon_s = ee.Number(epsilon_s)
epsilon_w = ee.Number(epsilon_w)

t_threshold = ee.Number(t_threshold)


# --------------------------------------------------
# IMPORT IMAGE COLLECTIONS
# --------------------------------------------------

# Landsat 5 TM
imgCol_L5_TOA = ee.ImageCollection('LANDSAT/LT05/C01/T1')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', max_cloud_cover))\
    .select(['B6'])

imgCol_L5_SR = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', max_cloud_cover))\
    .map(fun_mask_ls_sr)

imgCol_L5_SR = imgCol_L5_SR.map(fun_bands_l57)

# Landsat 7 ETM+
imgCol_L7_TOA = ee.ImageCollection('LANDSAT/LE07/C01/T1')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', max_cloud_cover))\
    .select(['B6_VCID_2'])

imgCol_L7_SR = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', max_cloud_cover))\
    .map(fun_mask_ls_sr)

imgCol_L7_SR = imgCol_L7_SR.map(fun_bands_l57)

# Landsat 8 OLI-TIRS
imgCol_L8_TOA = ee.ImageCollection('LANDSAT/LC08/C01/T1')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', max_cloud_cover))\
    .select(['B10'])

imgCol_L8_SR = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))\
    .filter(ee.Filter.lt('CLOUD_COVER_LAND', max_cloud_cover))\
    .map(fun_mask_ls_sr)

imgCol_L8_SR = imgCol_L8_SR.map(fun_bands_l8)

# NCEP/NCAR Water Vapor Product
imgCol_WV = ee.ImageCollection('NCEP_RE/surface_wv')\
    .filterBounds(select_roi)\
    .filter(ee.Filter.calendarRange(year_start,year_end,'year'))\
    .filter(ee.Filter.calendarRange(month_start,month_end,'month'))


# --------------------------------------------------
# CALCULATE
# --------------------------------------------------

# TOA (Radiance) and SR
imgCol_L5_TOA = imgCol_L5_TOA.map(fun_radcal)
imgCol_L7_TOA = imgCol_L7_TOA.map(fun_radcal)
imgCol_L8_TOA = imgCol_L8_TOA.map(fun_radcal)

imgCol_L5_SR = ee.ImageCollection(join_l.apply(imgCol_L5_SR, imgCol_L5_TOA, maxDiffFilter))
imgCol_L7_SR = ee.ImageCollection(join_l.apply(imgCol_L7_SR, imgCol_L7_TOA, maxDiffFilter))
imgCol_L8_SR = ee.ImageCollection(join_l.apply(imgCol_L8_SR, imgCol_L8_TOA, maxDiffFilter))

imgCol_L5_SR = imgCol_L5_SR.map(fun_l_addband)
imgCol_L7_SR = imgCol_L7_SR.map(fun_l_addband)
imgCol_L8_SR = imgCol_L8_SR.map(fun_l_addband)

# Water Vapor
imgCol_L5_SR = ee.ImageCollection(join_wv.apply(imgCol_L5_SR, imgCol_WV, maxDiffFilter))
imgCol_L7_SR = ee.ImageCollection(join_wv.apply(imgCol_L7_SR, imgCol_WV, maxDiffFilter))
imgCol_L8_SR = ee.ImageCollection(join_wv.apply(imgCol_L8_SR, imgCol_WV, maxDiffFilter))

imgCol_L5_SR = imgCol_L5_SR.map(fun_wv_scale)
imgCol_L7_SR = imgCol_L7_SR.map(fun_wv_scale)
imgCol_L8_SR = imgCol_L8_SR.map(fun_wv_scale)

# Atmospheric Functions
imgCol_L5_SR = imgCol_L5_SR.map(fun_af1(cs_l5))
imgCol_L5_SR = imgCol_L5_SR.map(fun_af2(cs_l5))
imgCol_L5_SR = imgCol_L5_SR.map(fun_af3(cs_l5))

imgCol_L7_SR = imgCol_L7_SR.map(fun_af1(cs_l7))
imgCol_L7_SR = imgCol_L7_SR.map(fun_af2(cs_l7))
imgCol_L7_SR = imgCol_L7_SR.map(fun_af3(cs_l7))

imgCol_L8_SR = imgCol_L8_SR.map(fun_af1(cs_l8))
imgCol_L8_SR = imgCol_L8_SR.map(fun_af2(cs_l8))
imgCol_L8_SR = imgCol_L8_SR.map(fun_af3(cs_l8))

# Delta and Gamma Functions
imgCol_L5_SR = imgCol_L5_SR.map(fun_delta_l5)
imgCol_L7_SR = imgCol_L7_SR.map(fun_delta_l7)
imgCol_L8_SR = imgCol_L8_SR.map(fun_delta_l8)

imgCol_L5_SR = imgCol_L5_SR.map(fun_gamma_l5)
imgCol_L7_SR = imgCol_L7_SR.map(fun_gamma_l7)
imgCol_L8_SR = imgCol_L8_SR.map(fun_gamma_l8)

# Merge Collections
imgCol_merge = imgCol_L8_SR.merge(imgCol_L7_SR).merge(imgCol_L5_SR)
imgCol_merge = imgCol_merge.sort('system:time_start')

# Parameters and Indices
imgCol_merge = imgCol_merge.map(fun_ndvi)
imgCol_merge = imgCol_merge.map(fun_ndwi)
imgCol_merge = imgCol_merge.map(fun_tcg)
imgCol_merge = imgCol_merge.map(fun_tcb)
imgCol_merge = imgCol_merge.map(fun_tcw)

imgCol_merge = imgCol_merge.map(fun_fvc)
imgCol_merge = imgCol_merge.map(fun_epsilon_scale)
imgCol_merge = imgCol_merge.map(fun_epsilon)


# LST
imgCol_merge = imgCol_merge.map(fun_lst)
imgCol_merge = imgCol_merge.map(fun_mask_lst)


# --------------------------------------------------
# SPECTRAL TEMPORAL METRICS
# --------------------------------------------------
# Iterate over parameters and metrics
for parameter in select_parameters:

    # Mosaic imgCollection
    pseudodate = ee.Date('1960-01-01')
    subCol = ee.ImageCollection(imgCol_merge.select(parameter))
    dates = fun_getdates(subCol)
    ini_date = ee.Date(dates.get(0))
    ini_merge = subCol.filterDate(ini_date, ini_date.advance(1, 'day'))
    ini_merge = ini_merge.select(parameter).mosaic().set('system:time_start', ini_date.millis())
    ini = ee.List([ini_merge])
    imgCol_mosaic = ee.ImageCollection(ee.List(dates.iterate(fun_mosaic, ini)))
    imgCol_mosaic = imgCol_mosaic.map(fun_timeband)

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
                temp = temp.multiply(0.2)
                temp = temp.multiply(1000).int16()
            else:
                temp = temp.multiply(10000).int16()

        # Export to Drive
        filename = parameter+'_'+roi_filename+'_GEE_'+str(year_start)+'-'+str(year_end)+'_'+\
                   str(month_start)+'-'+str(month_end)+'_'+metric
        out = ee.batch.Export.image.toDrive(image=temp, description=filename,
                                            scale=pixel_resolution,
                                            maxPixels=1e13,
                                            region=select_roi['coordinates'][0],
                                            crs=epsg)
        process = ee.batch.Task.start(out)


'''
 References:
 
 Jiménez‐Muñoz, J.C.; Sobrino, J.A. A generalized single-channel method for retrieving
 land surface temperature from remote sensing data. Journal of Geophysical Research:
 Atmospheres 2003, 108.
 
 Jimenez-Munoz, J.C.; Cristobal, J.; Sobrino, J.A.; Soria, G.; Ninyerola, M.; Pons, X.;
 Pons, X. Revision of the Single-Channel Algorithm for Land Surface Temperature
 Retrieval From Landsat Thermal-Infrared Data. IEEE Transactions on Geoscience and
 Remote Sensing 2009, 47, 339–349.
 
 Jiménez-Muñoz, J.C.; Sobrino, J.A.; Skoković, D.; Mattar, C.; Cristóbal, J. Land
 Surface Temperature Retrieval Methods From Landsat-8 Thermal Infrared Sensor Data.
 IEEE Geoscience and Remote Sensing Letters 2014, 11, 1840–1843.
 
 Sobrino, J.A.; Jimenez-Munoz, J.C.; Soria, G.; Romaguera, M.; Guanter, L.; Moreno, J.;
 Plaza, A.; Martinez, P. Land Surface Emissivity Retrieval From Different VNIR and TIR
 Sensors. IEEE Transactions on Geoscience and Remote Sensing 2008, 46, 316–327.
'''

# ====================================================================================================#
# END
# ====================================================================================================#
