"""Script to transform an image from .nc file to geotiff
Authors= Ton Haer and Andres Diaz
"""

import fileinput
import string
import os
import gdal
import ogr
import osr
import gdalnumeric
import gdalconst

def nc_to_gtif(nc_input, nc_variable, gtiff_output, crs, nodatavalue):
    """Will transform a NetCDF to a GeoTIFF.
    Arguments:
        nc_input {string} -- string of data path to NetCDF file.
        nc_variable {string} -- variable from NetCDF.
        gtiff_out {string} -- GeoTIFF output path.
        res {float} -- original resolution of NetCDF.
        crs {int} -- EPSG number
    Returns:
        string -- stating converted NetCDF file
    """
    # import NetCDF as dataset
    nc_ds = nc.Dataset(nc_input, 'r', format='NETCDF4')
    # extract variable of interest
    ds_var = nc_ds.variables[nc_variable]
    # extract the fill variable
    var_FillValue = ds_var._FillValue
    # set masked values
    ds_var.set_auto_maskandscale(False)
    # flipup because of netcdf data convention
    ds_var = np.flipud(ds_var)
    # transform variables in numpy array
    ds_var = np.array(ds_var)
    # set fill variables to NaN
    ds_var[ds_var == var_FillValue] = np.nan
    # obtain dimensions of the input dataset
    ds_dim = ds_var.shape
    lon_col = int(ds_dim[1])
    lat_row = int(ds_dim[0])
    dim_lon = nc_ds.variables['lon']
    dim_lat = nc_ds.variables['lat']
    min_lon = np.min(dim_lon[:])
    max_lon = np.max(dim_lon[:])
    min_lat = np.min(dim_lat[:])
    max_lat = np.max(dim_lat[:])
    # determine lon en lat resolution
    res_lon = (max_lon - min_lon) / (lon_col - 1)
    res_lat = (max_lat - min_lat) / (lat_row - 1)
    # create GeoTIFF with appropriate dimensions at the output path
    _raster = gdal.GetDriverByName('GTiff').Create(gtiff_output, lon_col, lat_row, 1, gdal.GDT_Float32)
    # transform the GeoTIFF to match the NetCDF
    _raster.SetGeoTransform((min_lon-(0.5*res_lon), res_lon, 0, max_lat+(0.5*res_lat), 0, -res_lat))
    # set the projection of the GeoTIFF to the NetCDF projection
    source = osr.SpatialReference()
    source.ImportFromEPSG(crs)
    _raster.SetProjection(source.ExportToWkt())
    # write the NetCDF value to the created GeoTIFF
    _raster.GetRasterBand(1).SetNoDataValue(nodatavalue)
    _raster.GetRasterBand(1).WriteArray(ds_var, 0, 0)
    # close file
    _raster = None
    return basename(nc_input)+" is converted to GeoTIFF"

nc_to_gtif(nc_input, nc_variable, gtiff_output, crs, nodatavalue)

   BASEMAP = gdal.Open("Previous_Stored.tif") #an image to copy all the projections etc..
   #band_INUN_Depths = BASEMAP.GetRasterBand(1) #only if we want to see the band f that tif
   #INUNDATION_Y = (band_INUN_Depths.ReadAsArray(0, 0, INUNDATION_Depths.RasterXSize, INUNDATION_Depths.RasterYSize).astype(np.float64))

#_______________to safe the 




    nc_file = nc.Dataset(nc_to_covert, 'r', format='NETCDF4');
    variable=nc_file.variables['name_of_the_variable']
    variable.set_auto_maskandscale(False)
    variable = np.flipud(cell_size_km2) #in the nc format is upside down

    outDs = gdal.GetDriverByName('GTiff').Create("Name_OF_output.tif", BASEMAP.RasterXSize, BASEMAP.RasterYSize, gdal.GDT_Float64)
    outDs.SetGeoTransform(BASEMAP.GetGeoTransform())  #we make the same Transform than the base image
    outDs.SetProjection(BASEMAP.GetProjection())
    outDs.GetRasterBand(1).WriteArray(variable,0,0)


#_____________________________________________________________if there is no a tiff twin to copy transformation, cells size tc..

# 1. Define pixel_size and NoData value of new raster
NoData_value = -9999
x_res = 1 # assuming these are the cell sizes
y_res = 1 # change as appropriate
pixel_size = 1

# 2. Filenames for in- and output
_in = r"C:\Users\HP\Desktop\hcmc\Flood_depths\6July2016_Elevation_measure_230cm\PointsExa.shp"
_out = r"C:\Users\HP\Desktop\hcmc\Flood_depths\6July2016_Elevation_measure_230cm\exa.tif"

# 3. Open Shapefile
source_ds = ogr.Open(_in)
source_layer = source_ds.GetLayer()
x_min, x_max, y_min, y_max = source_layer.GetExtent()

list_values=[]
list_x=[]
list_y=[]

for i in range (0,source_layer.GetFeatureCount()):
	FeaturePos=source_layer.GetFeature(i)
	X_Value=FeaturePos.GetField(0)
	Y_Value=FeaturePos.GetField(1)
	valuePos=FeaturePos.GetField(2)
	list_values.append(valuePos)
	list_x.append(X_Value)
	list_y.append(Y_Value)



# 4. Create Target - TIFF
cols = int( (x_max - x_min) / x_res )
rows = int( (y_max - y_min) / y_res )

_raster = gdal.GetDriverByName('GTiff').Create(_out, cols, rows, 1, gdal.GDT_Byte)
_raster.SetGeoTransform((x_min, x_res, 0, y_max, 0, -y_res))
_band = _raster.GetRasterBand(1)
_band.SetNoDataValue(NoData_value)

print lyr.GetMetadata()
# 5. Rasterize why is the burn value 0... isn't that the same as the background?
exa=gdal.RasterizeLayer(_raster, [1], source_layer, burn_values=[0])
exa=gdal.RasterizeLayer(_raster, [1], source_layer, options = ["ALL_TOUCHED=TRUE"])
exa=gdal.RasterizeLayer(_raster, [1], source_layer, options = ["ALL_TOUCHED=TRUE", "ATTRIBUTE=POLYNUMB"])

err = gdal.RasterizeLayer(rasterDS, [1], shpLayer, options = ["ALL_TOUCHED=TRUE", "ATTRIBUTE=POLYNUMB"])

