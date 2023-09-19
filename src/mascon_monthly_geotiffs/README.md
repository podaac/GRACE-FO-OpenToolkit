# GRACE

## Script
write_singlemonth_geotif.py

## Purpose:
This Python script is used to convert the JPL GRACE Mascon file from netCDF4 to GeoTIFF format. The most current Mascon dataset can be found and downloaded from https://podaac.jpl.nasa.gov/dataset/TELLUS_GRFO_L3_JPL_RL06.1_LND_v04?ids=&values=&search=rl06.1v04&provider=POCLOUD. This script decomposes the multi-year monthly Mascon netCDF file into single GeoTIFF files for each month.  

## Inputs and Outputs:
  1. Inputs: data path and filename (.nc)
  2. Outputs: output path.

Users need to manually modify the code at line 26 and 27 for the input and output strings. Output *filenames* are auto created. By default, the script will look in a local 'input' directory (./input/...) and write files to a local output (./output) directory.

update:

```
input_mascon_file = Path('./input/GRCTellus.JPL.200204_202008.GLO.RL06M.MSCNv02CRI.nc')
output_geotiff_path = Path('./output')
```
to

```
input_mascon_file = Path('PATH/TO/YOUR/INPUTFILE')
output_geotiff_path = Path('/PATH/TO/YOUR/OUTPUTFILE')
```

## Pre-requisites:
  The code works with both Python 2 and Python 3. The code has been tested on Mac OS and PC (Windows).
  Required Python libraries are listed below:
    1. netCDF - https://unidata.github.io/netcdf4-python/netCDF4/index.html
    2. gdal - https://pypi.org/project/GDAL/
    3. pathlib - https://www.python.org/dev/peps/pep-0428/

## PO.DAAC contacts:
*  munish.sikka@jpl.nasa.gov
*  wen-hao.li@jpl.nasa.gov
*  Help Desk: podaac@podaac.jpl.nasa.gov
