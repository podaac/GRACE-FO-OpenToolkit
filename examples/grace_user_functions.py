import requests
import s3fs
import json, copy
import numpy as np
from osgeo import gdal, ogr, osr
import time
from netCDF4 import Dataset, date2index
from datetime import date
from datetime import datetime as dt
import xarray as xr
import sys
import os.path
import subprocess
import netCDF4
import math
import csv
import fsspec
import pandas as pd
import statsmodels.api

def store_aws_keys(endpoint: str="https://archive.podaac.earthdata.nasa.gov/s3credentials"):    
    with requests.get(endpoint, "w") as r:
        accessKeyId, secretAccessKey, sessionToken, expiration = list(r.json().values())

    creds ={}
    creds['AccessKeyId'] = accessKeyId
    creds['SecretAccessKey'] = secretAccessKey
    creds['SessionToken'] = sessionToken
    creds['expiration'] = expiration
    
    return creds


def grace_connection(ShortName,grace_filename):
    #Source: Jinbo Wang (Email: jinbo.wang@jpl.nasa.gov)
    creds = store_aws_keys()
    #print(creds)
    s3 = s3fs.S3FileSystem(
    key = creds['AccessKeyId'],
    secret = creds['SecretAccessKey'],
    token = creds['SessionToken'],
    client_kwargs = {'region_name':'us-west-2'},
    )
    #print(f"\nThe current session token expires at {creds['expiration']}.\n")

# Ask PODAAC for the collection id using the 'short name'
    response = requests.get(
        url='https://cmr.earthdata.nasa.gov/search/collections.umm_json', 
        params={'provider': "POCLOUD",
                'ShortName': ShortName,
                'page_size': 1}
    )

    ummc = response.json()['items'][0]
    ccid = ummc['meta']['concept-id']
    #print(f'collection id: {ccid}')

    ss="podaac-ops-cumulus-protected/%s/*.nc"%ShortName
    GRACE_s3_files = np.sort(s3.glob(ss))
    full_filename=f'podaac-ops-cumulus-protected/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.1_V3/{grace_filename}'
    dataset = xr.open_dataset(s3.open(full_filename))
        
    return dataset

def read_grace_dataset(ShortName,grace_filename):
    dataset = grace_connection(ShortName,grace_filename)
    
    return dataset


def read_shapefile_singlelayer(shapefile,xdim,ydim):
    # Source: Jack McNelis (email: jack.mcnelis@jpl.nasa.gov)  
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(shapefile, 0)
    lyr = shp.GetLayer(0)
    ssrs = lyr.GetSpatialRef()
    wkt = ssrs.ExportToPrettyWkt()
    lyr.SetAttributeFilter("DN = 1")

    for i, feat in enumerate(lyr):
        if feat.GetField("DN") == '1':
            break

    feat = lyr.GetFeature(i)
        
    print(feat.GetField("DN"))  # Confirm the desired feature was selected before breaking the loop.
    geom = feat.GetGeometryRef()
    geojson = geom.ExportToJson()
    list(json.loads(geojson).keys())
    driver = ogr.GetDriverByName("MEMORY")
    featds = driver.CreateDataSource("MemoryDataset")
    newlyr = featds.CreateLayer("1", ssrs, geom_type=ogr.wkbPolygon)
    lyrid = ogr.FieldDefn("DN", ogr.OFTInteger)
    newlyr.CreateField(lyrid)
    lyrdefn = newlyr.GetLayerDefn()
    newfeat = ogr.Feature(lyrdefn)
    newgeom = ogr.CreateGeometryFromJson(geojson)
    newfeat.SetGeometry(newgeom)
    newfeat.SetField("DN", 1)
    newlyr.CreateFeature(newfeat)
    newfeat = None
    gt = (
    -180.0,                  # 0  X minimum (upper-left corner, the origin),
    360/xdim,                # 1  X resolution,
    0.0,                     # 2  X rotation,
    90.0,                    # 3  Y maximum (upper-left corner, the origin),
    0.0,                     # 4  Y rotation,
    -1*(180/ydim),           # 5  Y resolution
    )   
    print(gt)
    mask = gdal.GetDriverByName('MEM').Create(
    '',                   # No filename required for in-memory dataset.
    xdim, ydim,           # Match the dimensions of the GRACE global grid.
    1,                    # Output mask should contain only one band.
    gdal.GDT_Byte,        # Output type should be byte [0,1].
    )
   
    mask.SetGeoTransform(gt)      # Set the affine transform defined above as the mask's geotransform.

    mask.SetProjection(wkt)       # Set the wkt defn extracted from the shp as the target coordinate system.

    band = mask.GetRasterBand(1)  # Select the first and only band in raster mask.

    band.Fill(0)                  # Fill it with zeros.

    band.SetNoDataValue(0)        # Set its nodata value to zero.

    err = gdal.RasterizeLayer(
    mask,
    [1],                      # Set the target band(s); just the one band mask in this case.
    newlyr,                   # Set the source feature layer to rasterize in band 1.
    burn_values = [1],          # Fill the polygon coverage area with 1s.
    )

    mask.FlushCache()             # "Write" changes to the in-memory dataset.

    marr = mask.GetRasterBand(1).ReadAsArray()
    return marr
    
def read_shapefile_multilayers(shapefile,xdim,ydim,layer_name):
    # Source: Jack McNelis (email: jack.mcnelis@jpl.nasa.gov)     
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(shapefile, 0)
    lyr = shp.GetLayer()
    ssrs = lyr.GetSpatialRef()
    wkt = ssrs.ExportToPrettyWkt()
    for i, feat in enumerate(lyr):
        if feat.GetField("WMOBB_NAME") == layer_name:
            break
    feat = lyr.GetFeature(i)     
    #print(i)
    #print(feat.GetField("WMOBB_NAME"))
    #Get the feature's geojson representation:
    geom = feat.GetGeometryRef()
    geojson = geom.ExportToJson()
    list(json.loads(geojson).keys())
    driver = ogr.GetDriverByName("MEMORY")

    featds = driver.CreateDataSource("MemoryDataset")

    newlyr = featds.CreateLayer("COLORADO (also COLORADO RIVER)", ssrs, geom_type=ogr.wkbPolygon)

    lyrid = ogr.FieldDefn("ID", ogr.OFTInteger)

    newlyr.CreateField(lyrid)

    lyrdefn = newlyr.GetLayerDefn()

    newfeat = ogr.Feature(lyrdefn)

    newgeom = ogr.CreateGeometryFromJson(geojson)

    newfeat.SetGeometry(newgeom)

    newfeat.SetField("ID", 1)

    newlyr.CreateFeature(newfeat)

    newfeat = None 
    gt = (
    -180.0,                   # 0  X minimum (upper-left corner, the origin),
    360/xdim,                 # 1  X resolution,
    0.0,                      # 2  X rotation,
    90.0,                     # 3  Y maximum (upper-left corner, the origin),
    0.0,                      # 4  Y rotation,
    -1*(180/ydim),            # 5  Y resolution
    )
    
    #print(gt)
    mask = gdal.GetDriverByName('MEM').Create(
    '',                       # No filename required for in-memory dataset.
    xdim, ydim,               # Match the dimensions of the GRACE global grid.
    1,                        # Output mask should contain only one band.
    gdal.GDT_Byte,            # Output type should be byte [0,1].
    )
    
    mask.SetGeoTransform(gt)      # Set the affine transform defined above as the mask's geotransform.

    mask.SetProjection(wkt)       # Set the wkt defn extracted from the shp as the target coordinate system.

    band = mask.GetRasterBand(1)  # Select the first and only band in raster mask.

    band.Fill(0)                  # Fill it with zeros.

    band.SetNoDataValue(0)        # Set its nodata value to zero.

    err = gdal.RasterizeLayer(
    mask,
    [1],                      # Set the target band(s); just the one band mask in this case.
    newlyr,                   # Set the source feature layer to rasterize in band 1.
    burn_values = [1],          # Fill the polygon coverage area with 1s.
    )
    
    mask.FlushCache()             # "Write" changes to the in-memory dataset.

    marr = mask.GetRasterBand(1).ReadAsArray()
    
    # we need to return a bbox of shape
    env = feat.GetGeometryRef().GetEnvelope()
    bbox = [env[0], env[2], env[1], env[3]]
    return marr,bbox

#compute weighted area for this region_mask
def area(lats):
    # Modules:
    from pyproj import Geod
    # Define WGS84 as CRS:
    geod = Geod(ellps='WGS84')
    dx = 1/4.0 #mascon is half deg res so dx =1/4 
    c_area = lambda lat: geod.polygon_area_perimeter(np.r_[-dx,dx,dx,-dx], lat+np.r_[-dx,-dx,dx,dx])[0]
    out = []
    for lat in lats:
        out.append(c_area(lat))
    return np.array(out)
#source: https://github.com/podaac/the-coding-club/blob/main/notebooks/MEaSUREs-SSH-dask.ipynb

#first circshift/roll the mask for longitudes and then flip it to position 0-359 and S-N to match GRACE grid
#users can also adjust GRACE dataset to mask orientation, for this tutorial we keep datasets in GRACE grid format
def shift_to_GRACE_orientation(shift_lon,flip_lat,data_array,indexes_to_shift,axis_no):
    if shift_lon:
        temp_1a = np.roll(data_array, indexes_to_shift, axis =axis_no)
    else:
        temp_1a = copy.copy(data_array)
    if flip_lat:
        reoriented_grid = np.flipud(temp_1a)
    else:
        reoriented_grid = copy.copy(temp_1a)  
    return reoriented_grid

### Convert timemonth to year fraction
#Function toYearFraction(date) here converts time to year fraction. Since some of the datasets contain time as yearfrac, we convert all datasets time variable to consistent units.
def toYearFraction(date): #source: Internet:
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction
#https://stackoverflow.com/questions/6451655/how-to-convert-python-datetime-dates-to-decimal-float-years

def seamean_connection(ShortName):
    #Source: Jinbo Wang (Email: jinbo.wang@jpl.nasa.gov)
    creds = store_aws_keys()
#print(creds)
    s3 = s3fs.S3FileSystem(
    key = creds['AccessKeyId'],
    secret = creds['SecretAccessKey'],
    token = creds['SessionToken'],
    client_kwargs = {'region_name':'us-west-2'},
    )
    #print(f"\nThe current session token expires at {creds['expiration']}.\n")

# Ask PODAAC for the collection id using the 'short name'
    response = requests.get(
    url = 'https://cmr.earthdata.nasa.gov/search/collections.umm_json', 
    params = {'provider': "POCLOUD",
            'ShortName': ShortName,
            'page_size': 1}
    )

    ummc = response.json()['items'][0]
    ccid = ummc['meta']['concept-id']
    sea_ss = "podaac-ops-cumulus-protected/%s/*.txt"%ShortName
    sea_level_doclist = np.sort(s3.glob(sea_ss))
    sea_level_doc = sea_level_doclist[-1] #read latest file in the cloud for this dataset using last index of array
    sea_level_doc = "s3://" + sea_level_doc # reading txt file using fsspec and hence adding s3 bucket in path. 
    #print(sea_level_doc)
    with fsspec.open(sea_level_doc, mode="rb", anon=False, 
            key=creds["AccessKeyId"], secret=creds["SecretAccessKey"], 
            token=creds["SessionToken"]) as f:
        lines = f.readlines()
        all_data = [line.strip() for line in lines] #read all lines 
    return all_data
    

def read_sea_mean_doc(ShortName):
    all_data = seamean_connection(ShortName)
    cols_name_identifier = b'HDR column description'
    col_names = []
    linenumber_cols_desc_start = 0
    linenumber_cols_desc_end = 0
    linenumber_hdr_ends = 0
    for i in range(0,len(all_data)):
        current_line = all_data[i]
        if current_line==cols_name_identifier:
            linenumber_cols_desc_start = i+1;#next line contains first col desc
            break

    #first col desc is known. now determine total cols and end line for headers
    for i in range(linenumber_cols_desc_start,len(all_data)):
        current_line = all_data[i]
        if current_line==b'HDR':
            linenumber_cols_desc_end=i-1
            break
        
    for i in range(linenumber_cols_desc_end,len(all_data)):
        current_line = all_data[i]
        if current_line==b'HDR Header_End---------------------------------------':
            linenumber_hdr_ends = i
            break
    #print(linenumber_cols_desc_start)
    #print(linenumber_cols_desc_end)
    #print(linenumber_hdr_ends)
    data = all_data[linenumber_hdr_ends+1:]
    col_names = all_data[linenumber_cols_desc_start:linenumber_cols_desc_end+1]
    for col_counter in range(0,len(col_names)):
        if ('GMSL (Global Isostatic Adjustment (GIA) applied)' in     str(col_names[col_counter])):
            GMSL_col_no = col_counter
        if ('year+fraction' in str(col_names[col_counter])):
            yearfrac_col_no = col_counter
    #print(GMSL_col_no)
    #print(yearfrac_col_no)
    sea_level_data = np.loadtxt(all_data,skiprows=linenumber_hdr_ends+1)
    #print(sea_level_data.shape)
    sea_mean_vals = sea_level_data[:,GMSL_col_no]
    sea_time_vec = sea_level_data[:,yearfrac_col_no]
    return sea_mean_vals,sea_time_vec 

#Read the data, convert to desired units (mm) and convert time variable into yearfrac.
def read_thermosteric_data(ShortName,start_date, end_date):
    #Source: Jinbo Wang (Email: jinbo.wang@jpl.nasa.gov)
    creds = store_aws_keys()
#print(creds)
    s3 = s3fs.S3FileSystem(
    key = creds['AccessKeyId'],
    secret = creds['SecretAccessKey'],
    token = creds['SessionToken'],
    client_kwargs = {'region_name':'us-west-2'},
    )
    #print(f"\nThe current session token expires at {creds['expiration']}.\n")

# Ask PODAAC for the collection id using the 'short name'
    response = requests.get(
    url = 'https://cmr.earthdata.nasa.gov/search/collections.umm_json', 
    params = {'provider': "POCLOUD",
            'ShortName': ShortName,
            'page_size': 1}
    )

    ummc = response.json()['items'][0]
    ccid = ummc['meta']['concept-id']
    
    steric_ss = "podaac-ops-cumulus-protected/%s/*.nc"%ShortName
    thermosteric_s3_files = np.sort(s3.glob(steric_ss))
    steric_filename = thermosteric_s3_files[-1] #read latest file in the cloud for this dataset using last index of array
    #print(steric_filename)

    # reading thermosteric data 
    with xr.open_dataset(s3.open(steric_filename)) as steric_file_nc:
        monthly_steric_xr = steric_file_nc

    #extract this region from GRACE dataset from desired time period. 
    steric_data = monthly_steric_xr["thermosteric_ts"].sel(
    time = slice(start_date, end_date)).data
    steric_data = steric_data*1000 #convert to mm
    steric_time_xr = monthly_steric_xr["time"].sel(time=slice(start_date,end_date)).data
    steric_time = np.empty(steric_time_xr.shape[0])

    for month in range(steric_time_xr.shape[0]):
        steric_ts = pd.to_datetime(steric_time_xr[month])
        steric_time[month] = toYearFraction(steric_ts)
        
    return steric_data,steric_time

#User defined smooth function for non uniform x
def non_uniform_savgol(x, y, window, polynom):
    """
    Applies a Savitzky-Golay filter to y with non-uniform spacing
    as defined in x

    This is based on https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
    The borders are interpolated like scipy.signal.savgol_filter would do

    Parameters
    ----------
    x : array_like
        List of floats representing the x values of the data
    y : array_like
        List of floats representing the y values. Must have same length
        as x
    window : int (odd)
        Window length of datapoints. Must be odd and smaller than x
    polynom : int
        The order of polynom used. Must be smaller than the window size

    Returns
    -------
    np.array of float
        The smoothed y values
    """
    if len(x) != len(y):
        raise ValueError('"x" and "y" must be of the same size')

    if len(x) < window:
        raise ValueError('The data size must be larger than the window size')

    if type(window) is not int:
        raise TypeError('"window" must be an integer')

    if window % 2 == 0:
        raise ValueError('The "window" must be an odd integer')

    if type(polynom) is not int:
        raise TypeError('"polynom" must be an integer')

    if polynom >= window:
        raise ValueError('"polynom" must be less than "window"')

    half_window = window // 2
    polynom += 1

    # Initialize variables
    A = np.empty((window, polynom))     # Matrix
    tA = np.empty((polynom, window))    # Transposed matrix
    t = np.empty(window)                # Local x variables
    y_smoothed = np.full(len(y), np.nan)

    # Start smoothing
    for i in range(half_window, len(x) - half_window, 1):
        # Center a window of x values on x[i]
        for j in range(0, window, 1):
            t[j] = x[i + j - half_window] - x[i]

        # Create the initial matrix A and its transposed form tA
        for j in range(0, window, 1):
            r = 1.0
            for k in range(0, polynom, 1):
                A[j, k] = r
                tA[k, j] = r
                r *= t[j]

        # Multiply the two matrices
        tAA = np.matmul(tA, A)

        # Invert the product of the matrices
        tAA = np.linalg.inv(tAA)

        # Calculate the pseudoinverse of the design matrix
        coeffs = np.matmul(tAA, tA)

        # Calculate c0 which is also the y value for y[i]
        y_smoothed[i] = 0
        for j in range(0, window, 1):
            y_smoothed[i] += coeffs[0, j] * y[i + j - half_window]

        # If at the end or beginning, store all coefficients for the polynom
        if i == half_window:
            first_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    first_coeffs[k] += coeffs[k, j] * y[j]
        elif i == len(x) - half_window - 1:
            last_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    last_coeffs[k] += coeffs[k, j] * y[len(y) - window + j]

    # Interpolate the result at the left border
    for i in range(0, half_window, 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += first_coeffs[j] * x_i
            x_i *= x[i] - x[half_window]

    # Interpolate the result at the right border
    for i in range(len(x) - half_window, len(x), 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += last_coeffs[j] * x_i
            x_i *= x[i] - x[-half_window - 1]

    return y_smoothed
# Source: answer provided by user name 'Not a programmer' on this link 
# https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data