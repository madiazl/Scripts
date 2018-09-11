


# import packages
import os
import netCDF4 as nc
import numpy as np
import datetime
import gdal

os.chdir(os.getcwd())

def main(input_file, output_file,metadata):
    data = raster2array(input_file)
    fill_value = -9999 #getProps(input_file)
    #y = getProps(input_file)
    #x = getProps(input_file)
    prepare_nc(output_file, lat, lon, metadata)
    append_nc(output_file, variable, data[0:21600,0:43200], fill_value)


def raster2array(src_raster):
    file_cell_size_km2 = nc.Dataset(src_raster, 'r', format='h5');
    cell_size_km2=file_cell_size_km2.variables['data']
    cell_size_km2.set_auto_maskandscale(False)
    cell_size_km2 = np.flipud(cell_size_km2)
    return cell_size_km2

def getProps(src_raster):
    """
    This function retrieves properties from the geotiff input file
    Input:  raster geotiff file
    """
    # Open the file for reading
    raster = gdal.Open(src_raster)
    band = raster.GetRasterBand(1)
    # Retrieve properties
    fill_value  = band.GetNoDataValue()
    #y           = float(raster.RasterYSize)
    #x           = float(raster.RasterXSize)

    return fill_value

def prepare_nc(nc_file, lat, lon, metadata, units='Days since 2010-01-01 00:00:00', calendar='gregorian',Format="NETCDF4",zlib=True):
    """
    This function prepares a NetCDF file with given metadata, for a certain year
    The function assumes a gregorian calendar and a time unit 'Days since 2010-01-01 00:00:00'
    inputs:
        nc__file:     path to new netcdf file
        x:            xaxis
        y:            yaxis
        metadata:     dictionary with global attributes
        units:        time units to use in time axis
    """

    print('Setting up "' + nc_file + '"')
    start_time = datetime.datetime(2010, 1, 1)
    time       = nc.date2num(start_time, units=units, calendar=calendar)
    nc_trg     = nc.Dataset(nc_file, 'w', format=Format, zlib=True)

    print('Setting up dimensions and attributes. lat: ' + str(len(lat))+ ' lon: ' + str(len(lon)))
    nc_trg.createDimension('time', 1)
    nc_trg.createDimension('lat', len(lat))
    nc_trg.createDimension('lon', len(lon))

    time_var = nc_trg.createVariable('time','f4',('time',))
    time_var.units = units
    time_var.calendar = calendar
    time_var.standard_name = 'time'
    time_var.long_name = 'time'
    time_var.axis = 'T'
    time_var[:] = time

    y_var   = nc_trg.createVariable('lat','f4',('lat',))
    y_var.standard_name = 'latitude'
    y_var.long_name = 'latitude'
    y_var.units = 'degrees_north'
    y_var.axis = 'Y'

    x_var = nc_trg.createVariable('lon','f4',('lon',))
    x_var.standard_name = 'longitude'
    x_var.long_name = 'longitude'
    x_var.units = 'degrees_east'
    x_var.axis = 'X'
    
    y_var[:] = lat
    x_var[:] = lon

    #y_var[:] = np.arange(-90, -90 + y*1./120, 1./120)
    #x_var[:] = np.arange(-180, -180 + x*1./120, 1./120)

    projection= nc_trg.createVariable('projection','c')
    projection.long_name = 'wgs84'
    projection.EPSG_code = 'EPSG:4326'
    projection.proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    projection.grid_mapping_name = 'latitude_longitude'

    # add all attributes from user-defined metadata
    for attr in metadata:
        nc_trg.setncattr(attr, metadata[attr])
    nc_trg.sync()
    nc_trg.close()

    print('Finished preparing ' + nc_file)


def append_nc(nc_file, var_name, data, fill_value, dtype='f4', chunksizes=(1, 90, 180)):
    """
    Write a new variable to target NetCDF file.
    input:
        nc_file:        NetCDF object, referring to target file
        var_name:       String, variable name of source NetCDF
        data:           Data, belonging to the variable

    """

    # add the variable
    nc_obj = nc.Dataset(nc_file, 'a')
    variab = nc_obj.createVariable(var_name, dtype,
                                    ('time', 'lat', 'lon',),
                                    chunksizes=chunksizes,
                                    fill_value=fill_value,
                                    zlib=True)

    # add some general attributes usually used in lat lon data
    variab.units = 'POP/Cell'
    variab.long_name = 'Number of People per cell'
    variab.coordinates   = 'lat lon'

    # write data to the variable
    #data[np.isnan(data)] = fill_value
    variab[0, :, :] = data

    # close the file
    nc_obj.sync()
    nc_obj.close()




variable = 'Cell_Area'
print('Setup the lat/lon axes')
lon = np.arange(-180+1./240, 180, 1./120)
lat = np.arange(-90+1./240, 90, 1./120)

# List of global attributes in netCDF file
print('Making list of user-defined metadata')
metadata = {}
metadata['title'] = 'Population provided by WorldBank at 30 resolution'
metadata['institution']='WorldBank'
metadata['source'] = 'Jin Gao and Stephane Hallegatte'
metadata['resolution'] = '0.00833333 degrees x 0.00833333 degrees'
metadata['references'] = '-'
metadata['Conventions'] = '-'
metadata['comment'] = 'Data to cross the hazard maps of Aqueduct year 1'
metadata['summary'] = 'Worldwide POP'
metadata['history'] = 'File created: Andres Diaz'



outputfile='C:\Users\User1\Desktop\Lisa_Surfsara\cell_size_km2.nc'
inputfile='C:\Users\User1\Desktop\Lisa_Surfsara\cell_size_km2.h5'
main(inputfile, outputfile,metadata)