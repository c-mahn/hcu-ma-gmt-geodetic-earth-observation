import numpy as np
import math
import netCDF4 as nc
import os

legendreFactor1 = np.array([])
legendreFactor2 = np.array([])

def legendreFunctions(theta, maxDegree):
	global legendreFactor1
	global legendreFactor2
	
	if (legendreFactor1.shape[0] < maxDegree+1):
		print("Calculating Legendre Factors")
		legendreFactor1 = np.zeros([maxDegree+1,maxDegree+1],order='F')
		legendreFactor2 = np.zeros([maxDegree+1,maxDegree+1],order='F')
		legendreFactor1[1,1] = math.sqrt(3) 
		for n in range(2, maxDegree+1):
			legendreFactor1[n,n] = math.sqrt((2*n+1)/(2*n))
		for m in range(0, maxDegree):
			for n in range(m+1, maxDegree+1):
				f = (2*n+1)/((n+m)*(n-m))
				legendreFactor1[n,m] = math.sqrt(f*(2*n-1))
				legendreFactor2[n,m] = math.sqrt(f*(n-m-1)*(n+m-1)/(2*n-3))
				
				
	cosTheta = math.cos(theta)
	sinTheta = math.sin(theta)
	
	Pnm = np.zeros([maxDegree+1, maxDegree+1],order='F')
	Pnm[0,0] = 1
	for n in range(1, maxDegree+1):
		Pnm[n,n] = legendreFactor1[n,n]*sinTheta*Pnm[n-1,n-1]
	for m in range(0, maxDegree):
		Pnm[m+1,m] = legendreFactor1[m+1,m]*cosTheta*Pnm[m,m]
	for m in range(0, maxDegree):
		for n in range(m+2, maxDegree+1):
			Pnm[n,m] = legendreFactor1[n,m]*cosTheta*Pnm[n-1,m] - legendreFactor2[n,m]*Pnm[n-2,m]
			
	return Pnm


def save_global_grid(filename,grid):
   
    nlon = grid.shape[1]
    nlat = grid.shape[0]
    
    lon = np.arange(-math.floor(nlon/2),math.floor(nlon/2),1)
    lat = np.arange(-math.floor(nlat/2),math.floor(nlat/2),1)
    
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