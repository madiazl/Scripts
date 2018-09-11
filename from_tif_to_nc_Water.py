from osgeo import gdal
from osgeo import osr
import subprocess
import numpy as np
#from netcdf_funcs import *
import numpy as np
import netCDF4 as nc
from scipy.interpolate import interp1d
import pdb
from osgeo import gdal
import os
import csv
import pandas as pd
import subprocess
import datetime
import scipy.spatial as spatial
import matplotlib.pyplot as plt
from func_write_tif_nc import *
os.chdir(r'/projects/0/aqueduct/users/IVM/Geog_units')


gdal.AllRegister()

#read the filter water image

name_WFI=r"/projects/0/aqueduct/users/IVM/Geog_units/WATER_p50_1km_v2.tif"
file = gdal.Open(name_WFI)
band_file = file.GetRasterBand(1)
variable = (band_file.ReadAsArray(0, 0, file.RasterXSize, file.RasterYSize).astype(np.int16))
#we replace values the 65535 for NA values (-9999)
variabe=np.array(variable)
variable=np.flipud(variable)

x1 = np.linspace(-180+1./240, 180-1./240, 43200)
y1 = np.linspace(-90+1./240, 90-1./240, 21600)
time_list = [datetime.datetime(2016, 1, 1, 0, 0)]
RP='blank' #not used
metadata_global = {'history':'Timothy ask me to change the water raster from tof to nc', 'author': 'Andres Diaz, Philip Ward, Johanna Englhardt, Timothy Tiggeloven and grace of Good','source': 'GLOFRIS impact model','institution': 'Institute for Environmental Studies VU Amsterdam','title': 'Aqueduct Impact Model','references': 'http://floods.wri.org/','conventions': 'CF-1.6','project': 'Aqueduct Global Flood Analyzer'}
inun_index = ['Johanna!!!!!!!!!']
metadata_var = {'units': 'adimensional','standard_name': 'Water_pixels','long_name': 'Dirks raster water persistent','comment': 'Water mask elaborated for taking as a boundary filter'}


#________________________________________________________________________________________________________________________________________

filename=r"/projects/0/aqueduct/users/IVM/Geog_units/WATER_p50_1km_v2.nc"
write_nc = func_write_tif_nc(filename, x1, y1, time_list,RP, metadata_global,metadata_var,inun_index,variable)


