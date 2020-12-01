#!/usr/bin/env python3

# This python script reads the GRACE mascon Level 3 data available in netcdf
# 'GRCTellus.JPL.*.RL06M.MSCNv02CRI.nc' and creates geotiff for each time step.
#
# Each grace timestep has two geotiffs, one for variable lwe_thickness and other for uncertainty.
# Users need to download the GRACE mascon netcdf file and set its location in the variable input_mascon_file (line #25) ;
# Additionally users need to specify in line # 26 the location of output geotiff.
# For any comments and feedback, please contact: podaac@jpl.nasa.gov

import datetime
from pathlib import Path

import numpy as np

from netCDF4 import Dataset, num2date
from osgeo import gdal, osr

# Data can be downloaded from here: 
# https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06_V2

# You will need an EarthData account to access the PO.DAAC drive.

# replace the NetCDF (nc) filename with current mascon netcdf name as it gets updated every month
# input_mascon_file = Path('PATH/TO/YOUR/INPUTFILE')
# output_geotiff_path = Path('/PATH/TO/YOUR/OUTPUTFILE')

input_mascon_file = Path() / 'input' / 'GRCTellus.JPL.200204_202008.GLO.RL06M.MSCNv02CRI.nc' 
output_geotiff_path = Path() / 'output'
output_geotiff_path.mkdir(parents=True, exist_ok=True)

def gtiff_write(grid_data, data_res, output_geotiff_file, lat_bound_max, lon_bound_min):
    gdal_driver = gdal.GetDriverByName('GTiff')

    gtiff_name = str(output_geotiff_file)
    gtiff_out = gdal_driver.Create(gtiff_name, grid_data.shape[1], grid_data.shape[0], 1, gdal.GDT_Float32)

    out_geotransform = [lon_bound_min, data_res, 0, lat_bound_max, 0, -1*data_res]
    gtiff_out.SetGeoTransform(out_geotransform)
    
    out_wkt = osr.SpatialReference()
    out_wkt.ImportFromEPSG(4326)
    out_proj_wkt = out_wkt.ExportToWkt()
    gtiff_out.SetProjection(out_proj_wkt)
    
    gtiff_band_out = gtiff_out.GetRasterBand(1)
    gtiff_band_out.SetNoDataValue(-99999)
    gtiff_band_out.WriteArray(grid_data)

    gtiff_out = None
    gtiff_band_out = None


def gtiff_set_metadata(output_geotiff_file, ncdf, current_month_start_date, current_month_end_date, lon):
    attr_list = ncdf.ncattrs()
    attr_to_exclude = [
        'time_coverage_start',
        'time_coverage_end', 
        'geospatial_lon_min', 
        'geospatial_lon_max', 
        'date_created', 
        'months_missing'
    ]

    gtiff_name = str(output_geotiff_file)
    ras = gdal.Open(gtiff_name, gdal.GA_Update)

    for attr_name in attr_list:
        if attr_name not in attr_to_exclude: 
            ras.SetMetadataItem(attr_name, str(ncdf.getncattr(attr_name)))
    
    ras.SetMetadataItem('time_coverage_start', current_month_start_date)
    ras.SetMetadataItem('time_coverage_end', current_month_end_date)
    ras.SetMetadataItem('geospatial_lon_min', str(np.min(lon)))
    ras.SetMetadataItem('geospatial_lon_max', str(np.max(lon)))
    
    now = datetime.datetime.now()
    ras.SetMetadataItem('date_created', now.strftime("%Y-%m-%dT%H:%M:%S"))
    
    ras = None

def main():
    # this segment opens the netcdf file and reads its variables.
    ncdf = Dataset(input_mascon_file, mode='r')

    input_time = ncdf.variables['time'][:]
    lat = ncdf.variables['lat'][:]
    lon = ncdf.variables['lon'][:]
    lon_bounds= ncdf.variables['lon_bounds'][:]
    lat_bounds= ncdf.variables['lat_bounds'][:]
    input_lwe_thickness = ncdf.variables['lwe_thickness'][:] #time lat lon
    input_uncertainty = ncdf.variables['uncertainty'][:]

    time_bounds = ncdf.variables['time_bounds'][:]
    dtime = num2date(time_bounds[:], ncdf.variables['time_bounds'].units)

    # initialize the arrays of same size as input dataset
    grid_lwe = input_lwe_thickness * 0
    grid_uncertainty = input_uncertainty * 0
    timesize = len(input_time)

    data_res  = abs(lon[2 ] - lon[1]) #resolution of dataset
    indexes_to_shift=int(360 / (2 * data_res)) # no of global longitude divided by degree resolution; /2 gives mid point around which shift is done.

    #detremine if grid is south to north then flip the array else not 
    flip_lat = lat[0] < lat[-1]

    #detremine if longitudes starts at 0 and goes to n then circshift or roll the grid and longitude array  
    shift_lon = np.max(lon) > 180

    #this loop arrange dataset in 90 to -90 latitude and -180 to 180 longitude orientation and writes monthly geotiffs; 
    #works on array in the order: time*lat*lon ; 
    if input_lwe_thickness.shape == (len(input_time), len(lat), len(lon)):
        for time_index in range(0,timesize):
            temp_lwe = input_lwe_thickness[time_index,:,:]
            temp_uncertainty = input_uncertainty[time_index,:,:]
            
            # mascon netcdf contains grid : 0 to 359 longitudes and -90 to 90 degree latitudes. following lines sets the grid as 180 to -180 longitude and 90 to -90 latitude to be written into geotif 
            if flip_lat:
                temp_1a = np.flipud(temp_lwe)
                temp_1a_uncert = np.flipud(temp_uncertainty)
            else:
                temp_1a = temp_lwe
                temp_1a_uncert = temp_uncertainty
    
            if shift_lon:
                temp_1a = np.roll(temp_1a, indexes_to_shift, axis =1)
                temp_1a_uncert = np.roll(temp_1a_uncert, indexes_to_shift,axis =1)
                # arrange longitudes into -180 to 180 orientation 
                lon[np.where(lon>180)] = lon[np.where(lon > 180)] - indexes_to_shift
                lon = np.roll(lon, indexes_to_shift, axis=0)
                lon_bounds[np.where(lon_bounds >= 180)] = lon_bounds[np.where(lon_bounds >= 180)] - indexes_to_shift
                lon_bounds = np.roll(lon_bounds, indexes_to_shift, axis=0)
            
            grid_lwe[time_index,:,:] = temp_1a
            grid_uncertainty[time_index,:,:] = temp_1a_uncert
            lat_bound_max = np.max(lat_bounds)
            lon_bound_min = np.min(lon_bounds)
                    
            #determine the start and end period for each monthly timestep and use it in geotiff output filenames 
            start_year = dtime[time_index,0].year
            end_year = dtime[time_index,1].year
            
            start_day_of_year = dtime[time_index,0].timetuple().tm_yday
            end_day_of_year = dtime[time_index,1].timetuple().tm_yday
            
            start_timestring = str(start_year)  + ("%03d" %(start_day_of_year,))
            end_timestring = str(end_year)  + ("%03d" %(end_day_of_year,)) 
            
            current_month_start_date = dtime[time_index,0].strftime("%Y-%m-%dT%H:%M:%S")
            current_month_end_date = dtime[time_index,1].strftime("%Y-%m-%dT%H:%M:%S")
            
            #generate lwe_thickness tiff file 
            output_geotiff_file = 'mascon_lwe_thickness_' + start_timestring + '_' + end_timestring +  '.tif'
            output_geotiff_file = output_geotiff_path / output_geotiff_file

            # create geotiff
            gtiff_write(grid_lwe[time_index, :, :], data_res, output_geotiff_file, lat_bound_max, lon_bound_min)
            gtiff_set_metadata(output_geotiff_file, ncdf, current_month_start_date, current_month_end_date,lon)

            output_geotiff_uncert_file = 'mascon_uncertainty_' + start_timestring + '_' + end_timestring +  '.tif'
            output_geotiff_uncert_file = output_geotiff_path / output_geotiff_uncert_file

            # create geotiff
            gtiff_write(grid_uncertainty[time_index,:,:], data_res, output_geotiff_uncert_file, lat_bound_max, lon_bound_min)
            gtiff_set_metadata(output_geotiff_uncert_file, ncdf, current_month_start_date, current_month_end_date,lon)
    
    else:
        print('Array Dimensions Not in Desired Order; Time * Lat * Lon expected') 

    ncdf.close()    

if __name__ == '__main__':
    main()
