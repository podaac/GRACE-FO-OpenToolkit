# GRACE

## Script
write_singlemonth_geotif.py

## Purpose:
This Python script is used to convert the JPL GRACE Mascon file from netCDF4 to GeoTIFF format. The most current Mascon dataset can be found and downloaded from https://doi.org/10.5067/TEMSC-3JC62. This script decomposes the multi-year monthly Mascon netCDF file into single geoTiff file for each month.  

## Inputs and Outputs:
  1. Inputs: data path and filename (.nc)
  2. Outputs: output path.

Users need to manually modify the code at line 26 and 27 for the input and output string. Output filenames are auto created.

## Pre-requisites:
  The code works with both Python 2 and 3. Has been tested on Mac and PC.
  Required Python lib listed below:
    1. netCDF - https://unidata.github.io/netcdf4-python/netCDF4/index.html
    2. gdal - https://pypi.org/project/GDAL/ 
    3. pathlib - https://www.python.org/dev/peps/pep-0428/

## PO.DAAC contacts:
*  munish.sikka@jpl.nasa.gov
*  wen-hao.li@jpl.nasa.gov
*  Help Desk: podaac@podaac.jpl.nasa.gov
