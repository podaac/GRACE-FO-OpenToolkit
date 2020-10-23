# GRACE

## Script
write_singlemonth_geotif.py

## Purpose:
This Python script is used to convert the JPL GRACE Mascon file from netCDF4 to GeoTIFF format. The most current Mascon dataset can be found and downloaded from https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06_V2. It breaks the multi-year monthly Mascon netCDF file into single geoTiff file for each month.  

## Inputs and Outputs:
  1. Inputs: data path and filename (.nc)
  2. Outputs: output path.

User need to manually modify the code at line 26 and 27 for the input and output string. Output filenames are auto created.

## Pre-requisites:
  The code works with both Python 2 and 3. Has been tested on Mac and PC.
  Required Python lib listed below:
    1. netCDF
    2. gdal
    3. pathlib

## PO.DAAC contacts:
*  munish.sikka@jpl.nasa.gov
* wen-hao.li@jpl.nasa.gov
