# -*- coding: utf-8 -*-
"""
Created on  ...

@author: ...
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import functions


### Constants
# set constant for maximum degree
...
# set Earth gravity constant
...
# set Earth's radius
...

### read gravity field
print('read potential coefficients...')
# call your function to read the data of 'ITG-Grace2010s.gfc'
...
# limit the maximum degree of the potential coefficients to 40 (maxDegree)
...
...

# call your function to read the data of 'GRS80.gfc'
...
# subtract the potential coefficients of GRS80 from ITG-Grace2010s.gfc
...
...


### evaluate spherical harmonics on 1x1 grid
print('evaluate spherical harmonics on 1x1 grid...')

# define a vector with all longitudes from -180° to 180° (1-degree spacing)
...
# define a vetcor with all latitudes from -90° to 90° (1-degree spacing)
...
# convert latitude (phi) to co-latitude (theta) which is needed in the calculation
...

# initialize matrices for 
# 1) disturbing potential
# 2) gravity anomaly
# 3) gravity anomaly in satellite altitude
...
...
...
# define normal gravity
...

# define two radii:
# 1) for the Earth's surface
# 2) for satellite altitude
...
...

# In the following part, for each position on Earth (1-degree grid) the 
# disturbing potential, the gravity anomaly on Earth and the gravity 
# anomaly in satellite altitude is calculated.
# loop over all co-latitudes
for ...
    # calculate legendre functions for current theta
    ...
    # loop over all longitudes
    for ...
        # initialize spherical harmonic sum with zero
        ...
        # loop over all degrees
        for ...
            # calculate degree-dependent factors for
            # 1) disturbing potential
            # 2) gravity anomaly
            # 3) gravity anomaly in satellite altitude
            ...
            ...
            ...
            # loop over all orders from zero to current degree
            for ...
                # sum up spherical harmonics
                ...
            # apply degree-dependent factor and add sum to current value of
            # 1) disturbing potential
            # 2) gravity anomaly
            # 3) gravity anomaly in satellite altitude
            ...
            ...
            ...
            # reset spherical harmonic sum to zero for the next iteration
            ...

# multiply disturbing potential and geoid anomalies with leading factor
...
...
...
# convert disturbing potential to geoid heights
...

### display and save results
disp('display and save results...');

# you can complete the following command with your geoid heights
# or gravity anomalies to display the results in Python
plt.pcolor(N, cmap='RdBu_r') # <----- replace xxx with name of matrix to display
plt.colorbar()
plt.show()

# save disturbing potential and gravity anomalies to netCDF files for
# plotting with Generic Mapping Tools
functions.save_global_grid('geoid_height.nc',xxx1) # <----- replace xxx1 with name of matrix to save
functions.save_global_grid('grav_anom.nc',xxx2) # <----- replace xxx2 with name of matrix to save
functions.save_global_grid('grav_anom_sat.nc',xxx3) # <----- replace xxx3 with name of matrix to save