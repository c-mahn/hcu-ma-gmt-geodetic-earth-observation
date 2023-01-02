# Higher Geodesy Homework
# #############################################################################

# Authors:
# Silas Teske
# Joshua Wolf
# Christopher Mahn

# #############################################################################

# Import of Libraries
# -----------------------------------------------------------------------------

# import string as st
# import random as r
# import re
# from turtle import position
import matplotlib.pyplot as plt
# from scipy import interpolate
import numpy as np
# import math as m
# import sys
import os
# from scipy.fft import fft, fftfreq
# from scipy import signal
import functions as fu

# Settings
# -----------------------------------------------------------------------------

verbose = True

# Constants
# -----------------------------------------------------------------------------
gravity_constant    = 3.986005000e+14   # m^3 / (kg * s^2)
normal_gravity      = 9.80665           # m / s^2
radius              = 6.378137000e+06   # m
height              = radius+450000     # m
degree_n            = 40                
rho_grad            = 180/np.pi

# Functions
# -----------------------------------------------------------------------------

def import_data_ITG(input_filename):
    """
    This function is used to the spherical harmonic coefficients c and s.

    Args:
        input_filename (str): This specifies the name of the file, that will be
        imported.
    """
    # Opening and reading file from disk
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(verbose):
        print(f'[Info] Opening file "{input_filename}"', end="\r")
    with open(os.path.join("data", input_filename)) as file:
        if(verbose):
            print(f'[Info] Reading file "{input_filename}"', end="\r")
        data = file.readlines()
    if(verbose):
        print(f'[Info] Read file "{input_filename}" successfully')

    # Formating the data from disk into a two-dimentional list
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    data = data[0:874]
    count = 0
    for i, e in enumerate(data):
        if(i >= 13):
            data[i] = {"L": int(e[5:8]), "M": int(e[9:12]), "C": float(e[13:32]), "S": float(e[33:52]), "sigmaC": float(e[53:72]), "sigmaS": float(e[73:92])}
        else:
            data[i] = None
            count += 1
    if(verbose):
        print(f'[Info] Detected file header with {count} lines')
    
    while True:
        if(data[0]==None):
            data.pop(0)
        else:
            break

    # Return the loaded data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(data)      

def import_data_GRS80(input_filename):
    """
    This function is used to the spherical harmonic coefficients c and s.

    Args:
        input_filename (str): This specifies the name of the file, that will be
        imported.
    """
    # Opening and reading file from disk
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(verbose):
        print(f'[Info] Opening file "{input_filename}"', end="\r")
    with open(os.path.join("data", input_filename)) as file:
        if(verbose):
            print(f'[Info] Reading file "{input_filename}"', end="\r")
        data = file.readlines()
    if(verbose):
        print(f'[Info] Read file "{input_filename}" successfully')

    # Formating the data from disk into a two-dimentional list
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    count = 0
    for i, e in enumerate(data):
        if(i >= 13):
            data[i] = {"L": int(e[5:8]), "M": int(e[9:12]), "C": float(e[13:34]), "S": float(e[35:54]), "sigmaC": float(e[55:74]), "sigmaS": float(e[75:94])}
        else:
            data[i] = None
            count += 1
    if(verbose):
        print(f'[Info] Detected file header with {count} lines')
    
    while True:
        if(data[0]==None):
            data.pop(0)
        else:
            break

    # Return the loaded data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(data)

# Assemble Matrix from loaded data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def assemble_matrix(data, index):
    if(verbose):
        print(f'[Info] Assembling the matrix for {index}')
    n = 0
    m = 0
    for i in data:
        if i["L"] >= n:
            n = i["L"]
        if i["M"] >= m:
            m = i["M"]
    matrix = np.zeros((n+1,m+1))

    for i in data:
        matrix[i["L"]][i["M"]] = i[index]
    return(matrix)

# Subtract the potential coefficients of GRS80 from ITG-Grace2010s.gfc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def subtract_matrix(matrix, data, index):
    if(verbose):
        print(f'[Info] Subtracting the potential coefficients of GRS80')
    for i in data:
        matrix[int(i["L"]), int(i["M"])] -= i[index]
    return(matrix)

# Evaluating spherical harmonics on 1° x 1° grid
# Calculating:
#   1) disturbing potential 
#   2) gravity anomaly on the surface
#   3) gravity anomaly in satellite altitude
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def calculate_spherical_harmonics(Cmatrix, Smatrix, long, colat):
    if(verbose):
        print(f'[Info] Evaluating spherical harmonics on 1° x 1° grid')
    T = np.zeros((len(colat), len(long)))
    Delta_g_surface = np.zeros((len(colat), len(long)))
    Delta_g_satellite = np.zeros((len(colat), len(long)))
    P_NxM = []

    # loop over all co-latitudes
    for i in colat:
        # calculate legendre functions for current theta
        P_NxM = fu.legendreFunctions(i/rho_grad, degree_n)
        
        # loop over all longitudes
        for j in long:
            print(f' Co-Lat: {int(i+1):03d}, Long: {int(j+1):+04d}', end="\r")
            # initialize spherical harmonic sum with zero
            T_sum = 0
            Delta_g_surface_sum = 0
            Delta_g_satellite_sum = 0
            Spherical_Harmonics_sum = 0
            
            # loop over all degrees
            for k in range(degree_n):
                T_n = 1**(k+1)
                Delta_g_surface_n = ((k-1)/radius) * 1**(k+1)
                Delta_g_satellite_n = ((k-1)/height) * (radius/height)**(k+1)
                for l in range(degree_n):
                    Spherical_Harmonics_sum += Cmatrix[k,l] * P_NxM[k,l] * np.cos((l*j/rho_grad)) + Smatrix[k,l] * P_NxM[k,l] * np.sin((l*j/rho_grad))
                
                # apply degree-dependent factor and add sum to current value of
                # 1) disturbing potential
                # 2) gravity anomaly
                # 3) gravity anomaly in satellite altitude
                T_sum += T_n * Spherical_Harmonics_sum
                Delta_g_surface_sum += Delta_g_surface_n * Spherical_Harmonics_sum
                Delta_g_satellite_sum += Delta_g_satellite_n * Spherical_Harmonics_sum
                
                # reset spherical harmonic sum to zero for the next iteration
                Spherical_Harmonics_sum = 0
            
            T[np.invert(int(i))][int(j) + 180] = T_sum
            Delta_g_surface[np.invert(int(i))][int(j) + 180] = Delta_g_surface_sum
            Delta_g_satellite[np.invert(int(i))][int(j) + 180] = Delta_g_satellite_sum
    
    # multiply disturbing potential and geoid anomalies with leading factor
    T_geoid_anomalies = (gravity_constant/radius) * T
    gravity_anomalies_surface = ((gravity_constant/radius) * Delta_g_surface) * 1000
    gravity_anomalies_satellite = ((gravity_constant/radius) * Delta_g_satellite) * 1000
    
    # convert disturbing potential to geoid heights
    N = T_geoid_anomalies / (gravity_constant/radius**2)
    return(N, gravity_anomalies_surface, gravity_anomalies_satellite)

# Beginning of the Main Programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    data_ITG = import_data_ITG(f"ITG-Grace2010s.gfc")
    data_GRS80 = import_data_GRS80(f"GRS80.gfc")
    # print(data_ITG[860]["C"])
    # print(data_GRS80[4]["L"])

    C_NxM = assemble_matrix(data_ITG, "C")
    S_NxM = assemble_matrix(data_ITG, "S")
    # print(C_NxM)
    # print(S_NxM)

    C_NxM_norm = subtract_matrix(C_NxM, data_GRS80, "C")
    S_NxM_norm = S_NxM # If all S values in GRS80 are 0, we can skip the substraction
    # S_NxM_norm = subtract_matrix(S_NxM, data_GRS80, "S")

    # Defining a vector with all longitudes from -180° to 180° (1-degree spacing)
    longitudes_vector = np.array(np.linspace(-180, 179, 360))
    # print(longitudes_vector)

    # Defining a vector with all co-latitudes from 0° to 180° (1-degree spacing)
    colatitudes_vector = np.array(np.linspace(0, 179, 180))
    # print(colatitudes_vector)

    # Evaluating spherical harmonics
    N, gravity_anomalies_surface, gravity_anomalies_satellite = calculate_spherical_harmonics(C_NxM_norm, S_NxM_norm, longitudes_vector, colatitudes_vector)

    plt.pcolor(N, cmap='RdBu_r')
    plt.colorbar()
    plt.show()

    fu.save_global_grid(os.path.join("data","geoid_height.nc"), N)
    fu.save_global_grid(os.path.join("data","grav_anom_surface.nc"), gravity_anomalies_surface)
    fu.save_global_grid(os.path.join("data","grav_anom_satellite.nc"), gravity_anomalies_satellite)

    # os.system(f'gmt begin geoid_height png')