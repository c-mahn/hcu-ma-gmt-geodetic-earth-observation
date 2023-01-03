# Geodetic Earth Observation - Task 02
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

def import_gfc(filename, ignore_until="end_of_head"):
    """
    This function is used to import gfc-files.

    Args:
        filename (str): This specifies the name of the file, that will be imported.
        ignore_until (str, optional): This is the key until everthing in the header will be ignored. Defaults to "end_of_head".
    """
    print(f'[Info] Importing file "{filename}" as gfc-file')
    ignore = True
    data = []
    for line in open(os.path.join("data", filename)):  # Open file and read line by line
        if(ignore):
            # Ignore lines until the line with the specified string is reached
            line = line.split(" ")[0]
            if(line == ignore_until):
                ignore = False
        else:
            line = line.split(" ")
            line = list(filter(None, line))  # Remove empty strings
            if(line[0] == "gfc"):            
                # Read the data
                line = {"L": int(line[1]),
                        "M": int(line[2]),
                        "C": float(line[3]),
                        "S": float(line[4]),
                        "sigma_C": float(line[5]),
                        "sigma_S": float(line[6])}
                data.append(line)
    return(data)


def assemble_matrix(data, value_index, coord_indices = ["L", "M"]):
    """
    This function is used to assemble a matrix from the data.

    Args:
        data (list): This is the data, that will be used to assemble the matrix.
        value_index (int/str): This is the index of the value, that will be used to assemble the matrix.
        coord_indices (list, optional): This . Defaults to ["L", "M"].
    """
    print(f'[Info] Assembling the matrix with dimensions {coord_indices[0]} x {coord_indices[1]}', end="\r")
    # Get the dimensions of the matrix
    dimension_x = 0
    dimension_y = 0
    for line in data:
        if line[coord_indices[0]] > dimension_x:
            dimension_x = line[coord_indices[0]]
        if line[coord_indices[1]] > dimension_y:
            dimension_y = line[coord_indices[1]]
    print(f'[Info] Assembling the matrix with dimensions {dimension_x+1}x{dimension_y+1} ({coord_indices[0]} x {coord_indices[1]})')
    matrix = np.zeros((dimension_x+1,dimension_y+1))

    # Assemble the matrix
    for line in data:
        matrix[line[coord_indices[0]]][line[coord_indices[1]]] = line[value_index]
    return(matrix)


def substract_matrices(matrix1, matrix2):
    """
    This function is used to substract two matrices with different dimmentions.
    The smaller matrix will be padded with zeros.
    """
    print(f'[Info] Substracting two matrices', end="\r")
    shape_1 = [len(matrix1), len(matrix1[0])]
    shape_2 = [len(matrix2), len(matrix2[0])]
    shape = [max(shape_1[0], shape_2[0]), max(shape_1[1], shape_2[1])]
    matrix = np.zeros((shape[0], shape[1]))
    shape = [len(matrix), len(matrix[0])]
    print(f'[Info] Substracting two matrices from {shape_1[0]}x{shape_1[1]} and {shape_2[0]}x{shape_2[1]} to a matrix with dimensions {shape[0]}x{shape[1]}')

    # Add the values of the first matrix
    for i in range(shape_1[0]):
        for j in range(shape_1[1]):
            matrix[i][j] += matrix1[i][j]

    # Substract the values of the second matrix
    for i in range(shape_2[0]):
        for j in range(shape_2[1]):
            matrix[i][j] -= matrix2[i][j]

    return(matrix)
            
    


# Beginning of the Main Programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    # Importing data
    # - - - - - - - -
    
    # Importing the data from Grace's gravity field model
    data_ITG = import_gfc("ITG-Grace2010s.gfc")
    
    # Importing the data from the normal gravity field model
    data_grs80 = import_gfc("GRS80.gfc")
    
    # Assembling the matrices
    itg_c = assemble_matrix(data_ITG, "C")
    itg_s = assemble_matrix(data_ITG, "S")
    itg_sigma_c = assemble_matrix(data_ITG, "sigma_C")
    itg_sigma_s = assemble_matrix(data_ITG, "sigma_S")

    data_grs80_c = assemble_matrix(data_grs80, "C")
    data_grs80_s = assemble_matrix(data_grs80, "S")
    data_grs80_sigma_c = assemble_matrix(data_grs80, "sigma_C")
    data_grs80_sigma_s = assemble_matrix(data_grs80, "sigma_S")
    
    # Substracting the matrices
    data_norm_c = substract_matrices(itg_c, data_grs80_c)
    data_norm_s = substract_matrices(itg_s, data_grs80_s)  # This step is unnecessary, because the s-values are always zero

    # Defining a vector with all longitudes from -180° to 180° (1-degree spacing)
    longitudes_vector = np.array(np.linspace(-180, 180, 361))

    # Defining a vector with all co-latitudes from 0° to 180° (1-degree spacing)
    colatitudes_vector = np.array(np.linspace(0, 180, 181))

    # Calculating the spherical harmonics
    
    