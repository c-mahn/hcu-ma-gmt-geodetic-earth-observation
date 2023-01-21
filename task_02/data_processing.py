# Geodetic Earth Observation - Task 02
# #############################################################################

# Authors:
# Silas Teske
# Joshua Wolf
# Christopher Mahn

# #############################################################################

# Import of Libraries
# -----------------------------------------------------------------------------

import main
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
mass                = 5.9722*10**24     # kg
degree_n            = 96
order_m             = 20                
rho_grad            = 180/np.pi
rho_water           = 1000              # kg / m^3  
grid_spacing        = 0.5

# Functions
# -----------------------------------------------------------------------------

def import_gfc(filename):
    """
    This function is used to import gfc-files.

    Args:
        filename (str): This specifies the name of the file, that will be imported.
        ignore_until (str, optional): This is the key until everthing in the header will be ignored. Defaults to "end_of_head".
    """
    # print(f'[Info] Importing file "{filename}" as gfc-file')
    ignore_until="end_of_head"
    ignore = True
    data = []
    for line in open(filename):  # Open file and read line by line
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


def import_gfc_from_folder(path):
    """
    This function is used to import all gfc-files from a folder.

    Args:
        path (str): This is the path to the folder, that contains the gfc-files.
    """
    # print(f'[Info] Importing all gfc-files from folder "{path}"')
    gfc_datasets = []
    for filename in os.listdir(path):
        if filename.endswith(".gfc"):
            data = import_gfc(os.path.join(path, filename))
            filename = filename.split(".")[0]
            date = filename.split("_")[-1]
            gfc_datasets.append({"date": date, "data": data})
    return(gfc_datasets)


def import_csv(filename, delimiter=";"):
    """
    This function imports a csv-file. If possible it will also make numbers into floats

    Args:
        filename (str): This is the filename, where to import die csv from.
        delimiter (str): This is the delimiter seperating the columns
    """
    with open(filename, "r") as file:
        content = file.readlines()
    data = []
    for line in content:
        line = line.split(delimiter)
        for index, entry in enumerate(line):
            try:
                line[index] = float(entry)
            except:
                line[index] = entry.strip()
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
    # print(f'[Info] Assemble matrix {coord_indices[0]} x {coord_indices[1]}', end="\r")
    # Get the dimensions of the matrix
    dimension_x = 0
    dimension_y = 0
    for line in data:
        if line[coord_indices[0]] > dimension_x:
            dimension_x = line[coord_indices[0]]
        if line[coord_indices[1]] > dimension_y:
            dimension_y = line[coord_indices[1]]
    # print(f'[Info] Assemble matrix {dimension_x+1}x{dimension_y+1} ({coord_indices[0]} x {coord_indices[1]})')
    matrix = np.zeros((dimension_x+1,dimension_y+1))
    matrix = matrix.tolist()

    # Assemble the matrix
    for line in data:
        matrix[line[coord_indices[0]]][line[coord_indices[1]]] = line[value_index]
    return(matrix)


def matrix_math(matrix1, matrix2, operator="+"):
    """
    This function is used to do a math operation on two different-sized matrices. The difference will be padded.
    """
    # print(f'[Info] Combining two matrices', end="\r")
    shape_1 = [len(matrix1), len(matrix1[0])]
    shape_2 = [len(matrix2), len(matrix2[0])]
    shape = [max(shape_1[0], shape_2[0]), max(shape_1[1], shape_2[1])]
    matrix = np.zeros((shape[0], shape[1]))
    matrix = matrix.tolist()
    shape = [len(matrix), len(matrix[0])]
    # print(f'[Info] Combining two matrices with shapes {shape_1[0]}x{shape_1[1]} {operator} {shape_2[0]}x{shape_2[1]} = {shape[0]}x{shape[1]}')

    # Add the values of the first matrix
    for i in range(shape_1[0]):
        for j in range(shape_1[1]):
            matrix[i][j] += matrix1[i][j]

    # Apply the operator
    for i in range(shape_2[0]):
        for j in range(shape_2[1]):
            if(operator == "+"):
                matrix[i][j] += matrix2[i][j]
            elif(operator == "-"):
                matrix[i][j] -= matrix2[i][j]
            elif(operator == "*"):
                matrix[i][j] *= matrix2[i][j]
            elif(operator == "/"):
                matrix[i][j] /= matrix2[i][j]
    return(matrix)


def calc_EWH(lamda, theta, cnm, snm, M, R, rho, k):
    """
    This function is used to calculate the EWH.

    Args:
        lamda (float): vector containing longitudes (in radian) for each position where the EWH shall be evaluated.
        theta (float): vector containing co-latitudes (in radian) for each position where the EWH shall be evaluated.
        cnm (list): triangular matrix containing the spherical harmonic coefficients from which the EWH shall be computed
        snm (list): triangular matrix containing the spherical harmonic coefficients from which the EWH shall be computed
        M (int): mass of the earth in kg
        R (float): radius of the earth in m
        rho (float): density of water (in kg/m3)
        k (float): vector of Load Love numbers until same degree as given cnm and snm.
    """
    # Converting from laura's format to our format
    lamda = list(set(lamda))
    theta = list(set(theta))

    max_degree = len(cnm)-1
    ewh = []
    for i in range(len(lamda)):
        ewh.append([])
        for j in range(len(theta)):
            ewh[i].append(0)
    for index_colatitude, colatitude in enumerate(theta):
        print(f'[{int(100*index_colatitude/len(theta)):03d}%]', end='\r')
        P = fu.legendreFunctions(colatitude, max_degree)  # Calculate the legendre functions
        for index_longitude, longitude in enumerate(lamda):
            sum_outer = 0
            for n in range(max_degree):
                sum_inner = 0
                for m in range(n):
                    sum_inner += cnm[n][m]*(P[n][m]*np.cos(m*longitude)) + snm[n][m]*(P[n][m]*np.sin(m*longitude))
                    # sum_inner += cnm[n][m]*(P[n][m]*np.cos(colatitude)*np.cos(m*longitude)) + snm[n][m]*(P[n][m]*np.cos(colatitude)*np.sin(m*longitude))
                sum_outer += (2*n+1)/(1+k[n]) * sum_inner
            ewh[index_longitude][index_colatitude] = (M/rho*4*np.pi*R**2) * float(sum_outer)

    # Convert matrix-data to laura's vector format
    vector = []
    for line in ewh:
        for entry in line:
            vector.append(entry)
    vector = np.array(vector)

    return(vector)


def gaussian_filtering_factors(degree, filter_radius=200):
    """
    This function is used to calculate the gaussian filtering factors.

    Args:
        degree (int): Maximum degree for which the factors are calculated.
        filter_radius (int, optional): Filter radius in kilometers. Defaults to 200000.
    """
    filter_radius = filter_radius * 1000
    b = (np.log(2)) / (1-np.cos(filter_radius/radius))
    w_0 = 1
    w_1 = (1+np.e**(-2*b)) / (1-np.e**(-2*b)) - (1/b)
    w = [w_0, w_1]
    while(len(w) < degree):
        w.append((-(2*len(w)-1) / b) * w[-1] + w[-2])
    np.array(w)
    return(w)


def apply_gaussian_filtering(matrix, degree="auto", filter_radius=200):
    """
    This function is used to apply the gaussian filtering to the given matrix.

    Args:
        matrix (list): Matrix which shall be filtered.
        degree (int): Maximum degree for which the factors are calculated. Defaults to "auto".
        filter_radius (int, optional): Filter radius in kilometers. Defaults to 200000.
    """
    if(degree == "auto"):
        degree = np.shape(matrix)[0]
    w = gaussian_filtering_factors(degree, filter_radius)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] *= w[i]
    return(matrix)


def latlon_from_polygon(polygon, resolution):
    data = fu.getGridfromPolygon(np.array(polygon), resolution)
    data = data.tolist()
    longitudes = []
    latitudes = []
    for i in data:
        longitudes.append(i[0])
        latitudes.append(i[1])
    longitudes = list(set(longitudes))
    latitudes = list(set(latitudes))
    return(latitudes, longitudes)
    


# Beginning of the Main Programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    # Importing all datasets
    # - - - - - - - - - - - -

    for index, dataset in enumerate(main.datasets):
        print(f'[Info][{index+1}/{len(main.datasets)}] Importing dataset: {dataset["name"]}')
        if(dataset["type"] == "gfc"):
            if(dataset["is_folder"]):
                data = import_gfc_from_folder(os.path.join(main.folder_data, dataset["name"]))
            else:
                data = import_gfc(os.path.join(main.folder_data, dataset["name"]))
        elif(dataset["type"] == "csv"):
            data = import_csv(os.path.join(main.folder_data, dataset["name"]), delimiter=dataset["delimiter"])
        elif(dataset["type"] == "single_column"):
            content = import_csv(os.path.join(main.folder_data, dataset["name"]), delimiter="\n")
            data = []
            for entry in content:
                data.append(entry[0])
        else:
            print(f'[Error] The dataset {dataset["name"]} has an unknown type')
            continue
        main.datasets[index]["data"] = data
    print(f'[Info] Importing datasets finished')


    # Create dataset of the augmented gravity field model of each month
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Create new dataset for the output
    output_datasets = {"name": "grace_augmented", "data": []}

    # Load the static grace observation model and assemble it into matrices
    itsg_grace_2018 = main.select_dataset(main.datasets, "name", "ITSG-Grace2018s.gfc")["data"]
    itsg_grace_2018_matrix_c = assemble_matrix(itsg_grace_2018, value_index="C")
    itsg_grace_2018_matrix_s = assemble_matrix(itsg_grace_2018, value_index="S")
    itsg_grace_2018_matrix_sigma_c = assemble_matrix(itsg_grace_2018, value_index="sigma_C")
    itsg_grace_2018_matrix_sigma_s = assemble_matrix(itsg_grace_2018, value_index="sigma_S")
    del itsg_grace_2018

    length = len(main.select_dataset(main.datasets, "name", "ITSG-Grace")["data"])  # Just for the progress bar
    for index, itsg_grace in enumerate(main.select_dataset(main.datasets, "name", "ITSG-Grace")["data"]):
        print(f'[Info][{index+1}/{length}] Creating augmented grace-dataset ({itsg_grace["date"]})', end="\r")  # Progress bar
        itsg_grace_dataset = itsg_grace["data"]

        # Assemble monthly grace observations into matrices
        itsg_grace_matrix_c = assemble_matrix(itsg_grace_dataset, value_index="C")
        itsg_grace_matrix_s = assemble_matrix(itsg_grace_dataset, value_index="S")
        itsg_grace_matrix_sigma_c = assemble_matrix(itsg_grace_dataset, value_index="sigma_C")
        itsg_grace_matrix_sigma_s = assemble_matrix(itsg_grace_dataset, value_index="sigma_S")

        # - static grace model
        itsg_grace_matrix_c = matrix_math(itsg_grace_matrix_c, itsg_grace_2018_matrix_c, operator="-")
        itsg_grace_matrix_s = matrix_math(itsg_grace_matrix_s, itsg_grace_2018_matrix_s, operator="-")
        itsg_grace_matrix_sigma_c = matrix_math(itsg_grace_matrix_sigma_c, itsg_grace_2018_matrix_c, operator="-")
        itsg_grace_matrix_sigma_s = matrix_math(itsg_grace_matrix_sigma_s, itsg_grace_2018_matrix_s, operator="-")

        # Assemble monthly grace coefficients into matrices
        deg1_dataset = []
        for dataset in main.select_dataset(main.datasets, "name", "deg1")["data"]:
            if(dataset["date"] == itsg_grace["date"]):  # Selecting the correct dataset
                deg1_dataset = dataset["data"]
                break
        deg1_matrix_c = assemble_matrix(deg1_dataset, value_index="C")
        deg1_matrix_s = assemble_matrix(deg1_dataset, value_index="S")
        deg1_matrix_sigma_c = assemble_matrix(deg1_dataset, value_index="sigma_C")
        deg1_matrix_sigma_s = assemble_matrix(deg1_dataset, value_index="sigma_S")

        # + monthly grace coefficients
        itsg_grace_matrix_c = matrix_math(itsg_grace_matrix_c, deg1_matrix_c, operator="+")
        itsg_grace_matrix_s = matrix_math(itsg_grace_matrix_s, deg1_matrix_s, operator="+")
        itsg_grace_matrix_sigma_c = matrix_math(itsg_grace_matrix_sigma_c, deg1_matrix_sigma_c, operator="+")
        itsg_grace_matrix_sigma_s = matrix_math(itsg_grace_matrix_s, deg1_matrix_sigma_s, operator="+")

        # Append the monthly augmented grace dataset
        output_datasets["data"].append({"date": itsg_grace["date"], "data": {"C": itsg_grace_matrix_c,
                                                                             "S": itsg_grace_matrix_s,
                                                                             "sigma_C": itsg_grace_matrix_sigma_c,
                                                                             "sigma_S": itsg_grace_matrix_sigma_s}})

    # Append all the datasets to the main datasets
    main.datasets.append(output_datasets)
    print(f'[Info][Done] Creating dataset of the augmented gravity field models')

    # Remove the temporary variables
    del output_datasets
    del itsg_grace_2018_matrix_c, itsg_grace_2018_matrix_s, itsg_grace_2018_matrix_sigma_c, itsg_grace_2018_matrix_sigma_s
    del itsg_grace_matrix_c, itsg_grace_matrix_s, itsg_grace_matrix_sigma_c, itsg_grace_matrix_sigma_s
    del deg1_matrix_c, deg1_matrix_s, deg1_matrix_sigma_c, deg1_matrix_sigma_s


    # Calculate the equivalent water height
    # -------------------------------------

    # Loading the love numbers
    love_numbers = np.array(main.select_dataset(main.datasets, "name", "loadLoveNumbers_Gegout97.txt")["data"])
    
    # Latutude and longitude of the region of interest
    latitudes_vector_rad, longitudes_vector_rad = latlon_from_polygon(main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"], 0.5)
    colatitudes_vector_rad = []
    for latitude in latitudes_vector_rad:
        colatitudes_vector_rad.append(np.pi/2 - latitude)

    # Converting to laura's format
    colatitudes_vector_laura = []
    longitudes_vector_laura = []
    for colatitude in latitudes_vector_rad:
        for longitude in longitudes_vector_rad:
            colatitudes_vector_laura.append(colatitude)
            longitudes_vector_laura.append(longitude)
    colatitudes_vector_laura = np.array(colatitudes_vector_laura)
    longitudes_vector_laura = np.array(longitudes_vector_laura)

    # Uncomment the following section to calculate the equivalent water height for each month
    ''' # Create new dataset for the equivalent water height of each month
    dataset_ewh = {"name": "ewh_monthly", "data": []}

    length = len(main.select_dataset(main.datasets, "name", "grace_augmented")["data"])  # Just for the progress bar
    for index, dataset in enumerate(main.select_dataset(main.datasets, "name", "grace_augmented")["data"]):
        print(f'[Info][{index+1}/{length}] Calculating equivalent water height ({dataset["date"]})', end="\r")  # Progress bar

        # Get the corresponding dataset
        dataset_date = dataset["date"]
        data = dataset["data"]

        # Converting the spherical harmonic coefficients into numpy-matrices
        cnm_matrix = np.array(data["C"])
        snm_matrix= np.array(data["S"])

        # Calculating the unfiltered equivalent water heights (ewh)

        # ewh = calc_EWH(longitudes_vector_rad, colatitudes_vector_rad, cnm_matrix, snm_matrix, mass, radius, rho_water, love_numbers[0:cnm_matrix.shape[0]])
        ewh = fu.calc_EWH_fast(longitudes_vector_rad, colatitudes_vector_rad, cnm_matrix, snm_matrix, mass, radius, rho_water, love_numbers[0:cnm_matrix.shape[0]])
        dataset_ewh["data"].append({"date": dataset_date, "data": ewh})

    # Append the equivalent water height dataset
    main.datasets.append(dataset_ewh)
    print(f'[Info][Done] Creating dataset of the equivalent water height')

    # Remove the temporary variables
    del love_numbers, longitudes_vector, colatitudes_vector, cnm_matrix, snm_matrix, ewh '''
    
    # Calculate equivalent water height of april 2008
    print(f'[Info] Calculating equivalent water height (2008-04)', end="\r")
    # Get the corresponding datasets
    cnm = np.array(main.select_dataset(main.select_dataset(main.datasets, "name", "grace_augmented")["data"], "date", "2008-04")["data"]["C"])
    snm = np.array(main.select_dataset(main.select_dataset(main.datasets, "name", "grace_augmented")["data"], "date", "2008-04")["data"]["S"])
    love_numbers_vector = love_numbers[0:np.array(main.select_dataset(main.select_dataset(main.datasets, "name", "grace_augmented")["data"], "date", "2008-04")["data"]["C"]).shape[0]]
    # Execute one of the following functions
    # ewh = calc_EWH(longitudes_vector_laura, colatitudes_vector_laura, cnm, snm, mass, radius, rho_water, love_numbers_vector)
    ewh = fu.calc_EWH_fast(longitudes_vector_laura, colatitudes_vector_laura, cnm, snm, mass, radius, rho_water, love_numbers_vector)

    # Convert the equivalent water height from laura's format into a dataset
    ewh_data = np.zeros((len(latitudes_vector_rad), len(longitudes_vector_rad)))
    ewh_data = ewh_data.tolist()
    for index, value in enumerate(ewh):
        ewh_data[int(index/len(longitudes_vector_rad))][index%len(longitudes_vector_rad)] = value

    main.datasets.append({"name": "ewh_2008-04",
                          "data": ewh_data,
                          "type": "list_of_lists"})
    ewh_csv = []
    for i, line in enumerate(main.select_dataset(main.datasets, "name", "ewh_2008-04")["data"]):
        for j, value in enumerate(line):
            ewh_csv.append([i, j, value])
    np.savetxt(os.path.join(main.folder_result, "2008-04_ewh.csv"), ewh_csv, delimiter=";")
    del ewh_csv, i, j, line, value
    print(f'[Info][Done] Creating dataset of the equivalent water height (2008-04)')

    # Filtering of the spherical harmonic coefficients
    # ------------------------------------------------
    
    filter_radius = 300
    selected_date = "2008-04"

    # Create new dataset for the filtered spherical harmonic coefficients
    print(f'[Info] Filtering the spherical harmonic coefficients', end="\r")
    grace_single_filtered = {"name": "gaussian_filter_coefficients"}

    # Select the grace_augmented dataset for April 2008
    grace_single = main.select_dataset(main.select_dataset(main.datasets, "name", "grace_augmented")["data"], "date", selected_date)["data"]

    # Filter the grace_augmented dataset for April 2008
    grace_single_filtered["C"] = apply_gaussian_filtering(grace_single["C"], filter_radius)
    grace_single_filtered["S"] = apply_gaussian_filtering(grace_single["S"], filter_radius)
    grace_single_filtered["date"] = grace_single["date"]
    
    print(f'[Info][Done] Creating dataset of the filtered spherical harmonic coefficients')



    # Computation of a time series of region averages
    # -----------------------------------------------













    # Interpolation of missing GRACE months
    # -----------------------------------------------

    # grace_missing_months = fu.interp_missing_months()

    # Estimation of the linear mass trend
    # -----------------------------------------------





    ''' 
    def calculate_spherical_Harmonics(norm_C, norm_S, longitudes, colatitudes):
        print(f'[Info] Evaluating spherical harmonics')
        T = np.zeros((len(colatitudes), len(longitudes)))
        delta_g_surface = np.zeros((len(colatitudes), len(longitudes)))
        delta_g_surface = delta_g_surface.tolist()
        delta_g_satellite = np.zeros((len(colatitudes), len(longitudes)))
        delta_g_satellite = delta_g_satellite.tolist()
        P = []

        # loop over all co-latitudes
        for colat in colatitudes:
            # calculate legendre functions for current theta
            P = fu.legendreFunctions(colat/rho_grad, degree_n)
            
            # loop over all longitudes
            for long in longitudes:
                print(f' Co-Latitude: {int(colat+1):03d}, Longitude: {int(long+1):+04d}', end="\r")
                # initialize spherical harmonic sum with zero
                T_sum = 0
                delta_g_surface_sum = 0
                delta_g_satellite_sum = 0
                spherical_Harmonics_sum = 0
                
                # loop over all degrees
                for n in range(degree_n):
                    T_n = 1**(n+1)
                    delta_g_surface_n = ((n-1)/radius) * 1**(n+1)
                    delta_g_satellite_n = ((n-1)/height) * (radius/height)**(n+1)

                    # loop over all orders
                    for m in range(order_m):
                        spherical_Harmonics_sum += norm_C[n,m] * P[n,m] * np.cos((m*long/rho_grad)) + norm_S[n,m] * P[n,m] * np.sin((m*long/rho_grad))
                    
                    # apply degree-dependent factor and add sum to current value of
                    # 1) disturbing potential
                    # 2) gravity anomaly
                    # 3) gravity anomaly in satellite altitude
                    T_sum += T_n * spherical_Harmonics_sum
                    delta_g_surface_sum += delta_g_surface_n * spherical_Harmonics_sum
                    delta_g_satellite_sum += delta_g_satellite_n * spherical_Harmonics_sum
                    
                    # reset spherical harmonic sum to zero for the next iteration
                    spherical_Harmonics_sum = 0
                
                T[np.invert(int(colat))][int(long) + 180] = T_sum
                delta_g_surface[np.invert(int(colat))][int(long) + 180] = delta_g_surface_sum
                delta_g_satellite[np.invert(int(colat))][int(long) + 180] = delta_g_satellite_sum
        
        # multiply disturbing potential and geoid anomalies with leading factor
        T_geoid_anomalies = (gravity_constant/radius) * T
        gravity_anomalies_surface = ((gravity_constant/radius) * delta_g_surface) * 1000
        gravity_anomalies_satellite = ((gravity_constant/radius) * delta_g_satellite) * 1000
        
        # convert disturbing potential to geoid heights
        N = T_geoid_anomalies / (gravity_constant/radius**2)
        return(N, gravity_anomalies_surface, gravity_anomalies_satellite)
        


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
        data_norm_c = matrix_math(itg_c, data_grs80_c, operator="-")
        data_norm_s = matrix_math(itg_s, data_grs80_s, operator="-")  # This step is unnecessary, because the s-values are always zero

        # Defining a vector with all longitudes from -180° to 180° (1-degree spacing)
        longitudes_vector = np.array(np.linspace(-180, 180, 361))

        # Defining a vector with all co-latitudes from 0° to 180° (1-degree spacing)
        colatitudes_vector = np.array(np.linspace(0, 180, 181))
    

        # Calculating the gravity anomalies
        # ---------------------------------

        # Defining the degree and order of the spherical harmonics
        T = np.zeros((len(colatitudes_vector), len(longitudes_vector)))
        delta_g_surface = np.zeros((len(colatitudes_vector), len(longitudes_vector)))
        delta_g_satellite = np.zeros((len(colatitudes_vector), len(longitudes_vector)))
        P_NxM = []

        # loop over all co-latitudes
        for i in colatitudes_vector:
            # calculate legendre functions for current theta
            P_NxM = fu.legendreFunctions(i/rho_grad, degree_n)
            
            # loop over all longitudes
            for j in longitudes_vector:
                print(f' Co-Lat: {int(i+1):03d}, Long: {int(j+1):+04d}', end="\r")
                # initialize spherical harmonic sum with zero
                T_sum = 0
                delta_g_surface_sum = 0
                delta_g_satellite_sum = 0
                spherical_Harmonics_sum = 0
                
                # loop over all degrees
                for k in range(degree_n):
                    T_n = 1**(k+1)
                    delta_g_surface_n = ((k-1)/radius) * 1**(k+1)
                    delta_g_satellite_n = ((k-1)/height) * (radius/height)**(k+1)
                    for l in range(degree_n):
                        spherical_Harmonics_sum += data_norm_c[k,l] * P_NxM[k,l] * np.cos((l*j/rho_grad)) + data_norm_s[k,l] * P_NxM[k,l] * np.sin((l*j/rho_grad))
                    
                    # apply degree-dependent factor and add sum to current value of
                    # 1) disturbing potential
                    # 2) gravity anomaly
                    # 3) gravity anomaly in satellite altitude
                    T_sum += T_n * spherical_Harmonics_sum
                    delta_g_surface_sum += delta_g_surface_n * spherical_Harmonics_sum
                    delta_g_satellite_sum += delta_g_satellite_n * spherical_Harmonics_sum
                    
                    # reset spherical harmonic sum to zero for the next iteration
                    spherical_Harmonics_sum = 0
                
                T[np.invert(int(i))][int(j) + 180] = T_sum
                delta_g_surface[np.invert(int(i))][int(j) + 180] = delta_g_surface_sum
                delta_g_satellite[np.invert(int(i))][int(j) + 180] = delta_g_satellite_sum
        
        # multiply disturbing potential and geoid anomalies with leading factor
        T_geoid_anomalies = (gravity_constant/radius) * T
        gravity_anomalies_surface = ((gravity_constant/radius) * delta_g_surface) * 1000
        gravity_anomalies_satellite = ((gravity_constant/radius) * delta_g_satellite) * 1000
        
        # convert disturbing potential to geoid heights
        N = T_geoid_anomalies / (gravity_constant/radius**2)

        plt.pcolor(N, cmap='RdBu_r')
        plt.colorbar()
        plt.show()

        fu.save_global_grid(os.path.join("data","geoid_height.nc"), N)
        fu.save_global_grid(os.path.join("data","grav_anom_surface.nc"), gravity_anomalies_surface)
        fu.save_global_grid(os.path.join("data","grav_anom_satellite.nc"), gravity_anomalies_satellite)
    '''
