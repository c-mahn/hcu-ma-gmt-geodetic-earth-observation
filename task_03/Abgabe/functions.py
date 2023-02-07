# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 16:41:02 2021

@author: Jensen
"""

import os
import netCDF4 as nc

def save_grid(filename,grid,lon,lat):
    ###########################################################################
    # Author: Laura Jensen
    # Date: 06.01.2021
    #
    # 
    # Description: Saves grid as netCDF file
    # 
    # Input:
    # filename         filename (with path) of grid to save (e.g. 'mdt.nc')
    #                  
    # grid             2D matrix with the grid to be saved
    # lon, lat         vectors containing the longitudes and latitudes of the
    #                  grid coordinates
    #
    ###########################################################################
    nlon = grid.shape[1]
    nlat = grid.shape[0]
    
    if os.path.exists(filename):
        os.remove(filename)
    
    f = nc.Dataset(filename,'w', format='NETCDF4') #'w' stands for write
    f.createDimension('lon', nlon)
    f.createDimension('lat', nlat)
    
    longitude = f.createVariable('lon', 'f4', 'lon')
    latitude = f.createVariable('lat', 'f4', 'lat')  
    data = f.createVariable('data', 'f4', ('lat','lon'))
    
    longitude[:] = lon #The "[:]" at the end of the variable instance is necessary
    latitude[:] = lat
    data[:,:] = grid
    
    f.close()