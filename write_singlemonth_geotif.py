# This python script reads the GRACE mascon Level 3 data available in netcdf
# 'GRCTellus.JPL.*.RL06M.MSCNv02CRI.nc' and creates geotiff for each time step.
#
# Each grace timestep has two geotiffs, one for variable lwe_thickness and other for uncertainty.
# Users need to download the GRACE mascon netcdf file and set its location in the variable input_mascon_file (line #26) ;
# Additionally users need to specify in line # 27 the location of output geotiff.
# For any comments and feedback, please contact: podaac@jpl.nasa.gov

import datetime
from datetime import date
from pathlib import Path

from netCDF4 import Dataset

import numpy as np
import pandas as pd

from osgeo import gdal, osr

# Data can be downloaded from here: 
# https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06_V2

# You will need an EarthData account to access the PO.DAAC drive.

# replace the NetCDF (nc) filename with current mascon netcdf name as it gets updated every month
input_mascon_file = Path() / 'input' / 'GRCTellus.JPL.200204_202008.GLO.RL06M.MSCNv02CRI.nc' 
output_geotiff_path = Path() / 'output'
output_geotiff_path.mkdir(parents=True, exist_ok=True)


def GTiff_Write(grid_data, sr, output_geotiff_file):
    Gdal_Driver = gdal.GetDriverByName('GTiff')
    GTiff_Name = str(output_geotiff_file)
    GTiff_Out = Gdal_Driver.Create(GTiff_Name, grid_data.shape[1], grid_data.shape[0], 1, gdal.GDT_Float32)
    Out_GeoTransform = [-180, sr, 0, 90, 0, -1*sr]
    GTiff_Out.SetGeoTransform(Out_GeoTransform)
    Out_Wkt = osr.SpatialReference()
    Out_Wkt.ImportFromEPSG(4326)
    Out_Proj_Wkt = Out_Wkt.ExportToWkt()
    GTiff_Out.SetProjection(Out_Proj_Wkt)
    GTiff_Band_Out = GTiff_Out.GetRasterBand(1)
    GTiff_Band_Out.SetNoDataValue(-99999)
    GTiff_Band_Out.WriteArray(grid_data)
    GTiff_Out = None
    GTiff_Band_Out = None


def GTiff_SetMetadata(output_geotiff_file, fh, current_month_start_date, current_month_end_date):
    ras = gdal.Open(str(output_geotiff_file), gdal.GA_Update)
    ras.SetMetadataItem('Conventions', str(fh.Conventions))
    ras.SetMetadataItem('Metadata_Conventions', str(fh.Metadata_Conventions))
    ras.SetMetadataItem('standard_name_vocabulary', str(fh.standard_name_vocabulary))
    ras.SetMetadataItem('title', str(fh.title))
    ras.SetMetadataItem('summary', str(fh.summary))
    ras.SetMetadataItem('keywords', str(fh.keywords))
    ras.SetMetadataItem('keywords_vocabulary', str(fh.keywords_vocabulary))
    ras.SetMetadataItem('platform', str(fh.platform))
    ras.SetMetadataItem('institution', str(fh.institution))
    ras.SetMetadataItem('creator_name', str(fh.creator_name))
    ras.SetMetadataItem('creator_email', str(fh.creator_email))
    ras.SetMetadataItem('creator_url', str(fh.creator_url))
    ras.SetMetadataItem('creator_type', str(fh.creator_type))
    ras.SetMetadataItem('creator_institution', str(fh.creator_institution))
    ras.SetMetadataItem('publisher_name', str(fh.publisher_name))
    ras.SetMetadataItem('publisher_email', str(fh.publisher_email))
    ras.SetMetadataItem('publisher_url', str(fh.publisher_url))
    ras.SetMetadataItem('publisher_type', str(fh.publisher_type))
    ras.SetMetadataItem('publisher_institution', str(fh.publisher_institution))
    ras.SetMetadataItem('project', str(fh.project))
    ras.SetMetadataItem('program', str(fh.program))
    ras.SetMetadataItem('id', str(fh.id))
    ras.SetMetadataItem('naming_authority', str(fh.naming_authority))
    ras.SetMetadataItem('source', str(fh.source))
    ras.SetMetadataItem('processing_level', str(fh.processing_level))
    ras.SetMetadataItem('acknowledgement', str(fh.acknowledgement))
    ras.SetMetadataItem('license', str(fh.license))
    ras.SetMetadataItem('product_version', str(fh.product_version))
    ras.SetMetadataItem('time_epoch', str(fh.time_epoch))
    ras.SetMetadataItem('time_coverage_start', str(current_month_start_date))
    ras.SetMetadataItem('time_coverage_end', str(current_month_end_date))
    ras.SetMetadataItem('geospatial_lat_min', str(fh.geospatial_lat_min))
    ras.SetMetadataItem('geospatial_lat_max', str(fh.geospatial_lat_max))
    ras.SetMetadataItem('geospatial_lat_units', str(fh.geospatial_lat_units))
    ras.SetMetadataItem('geospatial_lat_resolution', str(fh.geospatial_lat_resolution))
    ras.SetMetadataItem('geospatial_lon_min', str(-179.5))
    ras.SetMetadataItem('geospatial_lon_max', str(179.5))
    ras.SetMetadataItem('geospatial_lon_units', str(fh.geospatial_lon_units))
    ras.SetMetadataItem('geospatial_lon_resolution', str(fh.geospatial_lon_resolution))
    ras.SetMetadataItem('time_mean_removed', str(fh.time_mean_removed))
    ras.SetMetadataItem('postprocess_1', str(fh.postprocess_1))
    ras.SetMetadataItem('postprocess_2', str(fh.postprocess_2))
    ras.SetMetadataItem('postprocess_3', str(fh.postprocess_3))
    ras.SetMetadataItem('GIA_removed', str(fh.GIA_removed))
    ras.SetMetadataItem('C_20_substitution', str(fh.C_20_substitution))
    ras.SetMetadataItem('C_30_substitution', str(fh.C_30_substitution))
    ras.SetMetadataItem('user_note_1', str(fh.user_note_1))
    ras.SetMetadataItem('user_note_2', str(fh.user_note_2))
    now = datetime.datetime.now()
    ras.SetMetadataItem('date_created', str(now.strftime("%Y-%m-%dT%H:%M:%S")))
    ras = None

# this segment opens the netcdf file and reads its variables.

fh = Dataset(input_mascon_file, mode='r')

input_time = fh.variables['time'][:]
input_lwe_thickness = fh.variables['lwe_thickness'][:] #time lat lon
input_uncertainty = fh.variables['uncertainty'][:]

time_vals = fh.variables['time'][:]
time_bounds = fh.variables['time_bounds'][:]
time_epoch = date.toordinal(date(2002, 1, 1)) + 366

# initialize the arrays of same size as input dataset
grid_lwe = input_lwe_thickness * 0
grid_uncertainty = input_uncertainty * 0
timesize = len(input_time)

# this loop arrange dataset in 90 to -90 latitude and -180 to 180 longitude orientation and writes monthly geotiffs
for time_index in range(0, timesize):
    temp_lwe = input_lwe_thickness[time_index, :, :]
    temp_uncertainty = input_uncertainty[time_index, :, :]

    # mascon netcdf contains grid : 0 to 359 longitudes and -90 to 90 degree latitudes.
    # following lines sets the grid as 180 to -180 longitude and 90 to -90 latitude to be written into geotif
    temp_1a = np.flipud(temp_lwe)
    temp_1a_uncert = np.flipud(temp_uncertainty)
    
    grid_lwe[time_index, :, 360:720] = temp_1a[:, 0:360]
    grid_lwe[time_index, :, 0:360] = temp_1a[:, 360:720]
    
    grid_uncertainty[time_index, :, 360:720] = temp_1a_uncert[:, 0:360]
    grid_uncertainty[time_index, :, 0:360] = temp_1a_uncert[:, 360:720]

    # determine the start and end period for each monthly timestep and use it in geotiff output filenames
    start_time = time_bounds[time_index, 0] + time_epoch
    end_time = time_bounds[time_index, 1] + time_epoch
    
    start_timestamp = pd.to_datetime(start_time - 719529, unit='D')
    end_timestamp = pd.to_datetime(end_time - 719529, unit='D')
    
    start_year = start_timestamp.year
    end_year = end_timestamp.year
    
    start_day_of_year = start_timestamp.dayofyear
    end_day_of_year = end_timestamp.dayofyear
    
    start_timestring = str(start_year)  + ("%03d" %(start_day_of_year,))
    end_timestring = str(end_year)  + ("%03d" %(end_day_of_year,))
    
    current_month_start_date = str(start_timestamp.strftime("%Y-%m-%dT%H:%M:%S"))
    current_month_end_date = str(end_timestamp.strftime("%Y-%m-%dT%H:%M:%S"))
    
    # generate lwe_thickness tiff file
    
    # resolution of dataset is 0.5 degree
    sr = 0.5 

    output_geotiff_file = 'mascon_lwe_thickness_' + start_timestring + '_' + end_timestring +  '.tif'
    output_geotiff_file = output_geotiff_path / output_geotiff_file

    # create geotiff
    GTiff_Write(grid_lwe[time_index, :, :], sr, output_geotiff_file)
    GTiff_SetMetadata(output_geotiff_file, fh, current_month_start_date, current_month_end_date)

    output_geotiff_uncert_file = 'mascon_uncertainty_' + start_timestring + '_' + end_timestring +  '.tif'
    output_geotiff_uncert_file = output_geotiff_path / output_geotiff_uncert_file

    # create geotiff
    GTiff_Write(grid_uncertainty[time_index,:,:], sr, output_geotiff_uncert_file)
    GTiff_SetMetadata(output_geotiff_uncert_file, fh, current_month_start_date, current_month_end_date)

fh.close()
