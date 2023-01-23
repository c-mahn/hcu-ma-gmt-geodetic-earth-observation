import numpy as np
import math
import netCDF4 as nc
import os
from matplotlib import path

legendreFactor1 = np.array([])
legendreFactor2 = np.array([])

def legendreFunctions(theta, maxDegree):
	global legendreFactor1
	global legendreFactor2
	
	if (legendreFactor1.shape[0] < maxDegree+1):
		# print("[Info] Calculating Legendre Factors")
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

def readData(filename,numskip):
    data = np.loadtxt(filename,skiprows=numskip,usecols=(1,2,3,4))
    maxDeg = int(np.max(data[:,0]))
    cnm = np.zeros([maxDeg+1,maxDeg+1],order='F')
    snm = np.zeros([maxDeg+1,maxDeg+1],order='F')

    imax = np.shape(data)[0]
    for i in range(imax):        
        cnm[int(data[i,0]),int(data[i,1])] = data[i,2]
        snm[int(data[i,0]),int(data[i,1])] = data[i,3]

    return cnm, snm


def getGridfromPolygon(poly,gridSpacing):
    ###########################################################################
    # Author: Laura Jensen
    # Date: 01.12.2020
    # last change: 02.12.2020
    # 
    # Description: Returns grid cell coordinates and area weights for all
    # points inside the given polygon
    # 
    # Input:
    # poly             matrix with two columns, the longitudes and latitudes
    #                  (in degree) of the corner coordinates of the polygon 
    # grid_spacing     grid spacing (in degree) of the resulting grid
    #                  coordinates
    # 
    # Output:
    # grid_region_area   matrix with three columns, the longitudes and
    #                    latitudes (in degree) of all points within the
    #                    given polygon and the area weight of the grid cell
    ###########################################################################
    
    poly_path = path.Path(poly[:,0:2])
    
    # define vectors containing all coordinates from a global grid
    lon = np.arange(-180+gridSpacing/2,180+gridSpacing/2,gridSpacing)
    lat = np.arange(-90+gridSpacing/2,90+gridSpacing/2,gridSpacing)
    lonM = np.tile(lon[:,np.newaxis],[1,len(lat)])
    latM = np.tile(lat,[len(lon),1])
    lonV = lonM.flatten()
    latV = latM.flatten()
    
    grid_glo = np.hstack((lonV[:,np.newaxis], latV[:,np.newaxis]))
    
    # find out which points are in the polygon
    in_poly = poly_path.contains_points(grid_glo)
    # get the indices of the points in the polygon
    idx = np.where(in_poly)
    # extract the coordinates of the global grid within the polygon
    grid_region = grid_glo[idx]
    
    # get area weights
    deltaLamda = gridSpacing*math.pi/180.0
    deltaPhi = gridSpacing*math.pi/180.0
    phi = grid_region[:,1]*math.pi/180.0
    area = deltaLamda * 2.0 * np.sin(deltaPhi/2.0) * np.cos(phi)
    
    # stack coordinates and area together
    grid_region_area = np.hstack((grid_region,area[:,np.newaxis]))
    
    # sort coordinates by latitude because this accelerates runtime in
    # calc_EWH_fast
    grid_region_area = grid_region_area[grid_region_area[:,1].argsort(),:]
    
    return grid_region_area


def calc_EWH_fast(lamda,theta,cnm,snm,M,R,rho,k):
    ###########################################################################
    # Author: Laura Jensen
    # Date: 28.10.2018
    # last change: 01.12.2020
    # 
    # Description: Calculates equivalent water heights (EWH) from given
    # spherical harmonic coefficients of the gravitational potential for every
    # given position
    # 
    # Input:
    # lamda, theta     vectors containing the longitudes and 
    #                  co-latitudes (in radian) for each position where the
    #                  EWH shall be evaluated. 
    # cnm, snm         triangular matrix containing the spherical harmonic 
    #                  coefficients from which the EWH shall be computed 
    # M, R             mass (in kg) and radius (in m) of the Earth
    # rho			   density of water (in kg/m^3)
    # k                vector of load love numbers until 
    #                  same degree as given cnm and snm
    # 
    # Output:
    # ewh              vector with same length as lambda
    #                  and theta containing the EWH
    ###########################################################################
    
    # define maximum degree
    maxDegree = cnm.shape[0]-1
    # initialize EWH
    ewh = np.zeros(len(theta))
    
    # calculate degree dependent factor
    degree = np.arange(0,maxDegree+1,1)
    radialFactor = (2*degree+1)/(1 + k)
    # calculate leading factor
    factor = M/(rho*4*math.pi*(R**2))
    
    # order is needed in loop
    order = np.arange(0,maxDegree+1,1)
    theta_alt = None
    
    # loop over all positions
    for i in range(len(theta)):
        Cnm = np.zeros((maxDegree+1,maxDegree+1))
        Snm = np.zeros((maxDegree+1,maxDegree+1))
        cosLambda = np.cos(order*lamda[i])
        sinLambda = np.sin(order*lamda[i])
        # only if theta changes, calc new legendre functions
        if (theta[i] != theta_alt):
            Pnm = legendreFunctions(theta[i],maxDegree)

        # multiply columns of legendre functions with cos/sin
        for j in range(Pnm.shape[1]):
            Cnm[:,j] = Pnm[:,j]*cosLambda[j]
            Snm[:,j] = Pnm[:,j]*sinLambda[j]

        theta_alt = theta[i]
        # sum over all orders (columns)
        Yn = np.sum(cnm*Cnm + snm*Snm,1)
        # multiply with degree dependent factor and sum over all degrees,
        # and multiply with leading factor
        ewh[i] = factor*np.sum(radialFactor*Yn)
    
    return ewh

def interp_missing_months(values):
    ###########################################################################
    # Author: Laura Jensen
    # Date: 28.10.2018
    # last change: 01.12.2020
    # 
    # Description: Interpolates the missing GRACE months in a given time
    # series (from 2003-01 to 2016-12) by taking the average of the previous 
    # and the successive month.
    # CAUTION: this function needs a time series starting at 2003-01 and ending
    # at 2016-12, otherwise the interpolation is not correct!!
    # 
    # Input:
    # values           column vector of time series with missing months
    # 
    # Output:
    # dec_yr           column vector of times for all months from 2003-01 to 
    #                  2016-12 in decimal years (without gaps)
    # values_interp    column vector of interpolated time series (without gaps)
    #                  with same length than dec_yr
    ###########################################################################

    # create date vector
    mn = []
    yr = []
    startyr = 2003;
    for i in range(14):
        mn = np.append(mn,np.arange(1,13))
        yr = np.append(yr,(startyr+i)*np.ones(12))
        
    # convert to decimal years
    dec_yr = yr+(mn/12)
    # missing months
    yr_miss = np.array([2003, 2011, 2011, 2012, 2012, 2013, 2013, 2013, 2014,
        2014, 2014, 2015, 2015, 2015, 2015, 2016, 2016, 2016])
    mn_miss = np.array([6, 1, 6, 5, 10, 3, 8, 9, 2, 7, 12, 5, 6, 10, 11, 4, 9, 10])

    # interpolation
    t_all = len(yr)
    values_interp = np.zeros(t_all)
    k = 0
    for i in range(t_all):
        # print(i)
        if np.any((yr[i]==yr_miss) * (mn[i]==mn_miss)):
            # print('interpolated')
            mean_val = np.mean(values[k:k+1])
            values_interp[i] = mean_val
        else:
            values_interp[i] = values[k]
            k = k+1
            
    return dec_yr, values_interp

