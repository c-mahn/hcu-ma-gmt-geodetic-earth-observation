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
import json

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
    data = np.array(data)
    return(data)


def import_data(dataset):
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
            data.append(float(entry[0]))
    else:
        print(f'[Error] The dataset {dataset["name"]} has an unknown type')
        return(None)
    return(data)


def export_data(dataset):
    try:
        if(type(dataset["data"]) is np.ndarray):
            try:
                longitudes = dataset["axis"][0]
                colatitudes = dataset["axis"][1]
                if(colatitudes[0] - colatitudes[1] == longitudes[0] - longitudes[1]):
                    spacing = colatitudes[1] - colatitudes[0]
                matrix = np.zeros((360, 180))
                for longitude in range(-180, 180, spacing):
                    for colatitude in range(0, 180, spacing):
                        for long_index, long in enumerate(longitudes):
                            for colat_index, col in enumerate(colatitudes):
                                if(int(long) == longitude and int(col) == colatitude):
                                    matrix[longitude][colatitude] = dataset["data"][long_index][colat_index]
                fu.save_global_grid(matrix, os.path.join("output", f'{dataset["name"]}.nc'))
            except:
                try:
                    with open(os.path.join("output", f'{dataset["name"]}.json'), 'w') as f:
                        json.dump(dataset, f)
                except:
                    os.remove(os.path.join("output", f'{dataset["name"]}.json'))
        elif(type(dataset["data"]) is list):
            for subdataset in dataset["data"]:
                export_data(subdataset)
    except:
        pass


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
    matrix = np.array(matrix)
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
        theta (float): vector containing co-longitudes (in radian) for each position where the EWH shall be evaluated.
        cnm (list): triangular matrix containing the spherical harmonic coefficients from which the EWH shall be computed
        snm (list): triangular matrix containing the spherical harmonic coefficients from which the EWH shall be computed
        M (int): mass of the earth in kg
        R (float): radius of the earth in m
        rho (float): density of water (in kg/m3)
        k (float): vector of Load Love numbers until same degree as given cnm and snm.
    """
    result = []
    max_degree = np.shape(cnm)[0]-1

    for i in range(len(lamda)):
        lamda_i = float(lamda[i])
        theta_i = float(theta[i])
        P = fu.legendreFunctions(theta_i, max_degree)  # Calculate the legendre functions
        sum_outer = 0
        for n in range(max_degree):
            sum_inner = 0
            for m in range(n):
                sum_inner += cnm[n][m]*(P[n][m]*np.cos(m*lamda_i)) + snm[n][m]*(P[n][m]*np.sin(m*lamda_i))
            sum_outer += (2*n+1)/(1+k[n]) * sum_inner
        result.append((M/rho*4*np.pi*R**2) * float(sum_outer))

    result = np.array(result)
    return(result)


def apply_ewh(dataset, M, R, rho, k, spacing=1, area=None):
    """
    This function is used to apply the EWH to a given dataset.
    
    Args:
        dataset (dict): Dataset containing the spherical harmonic coefficients.
        M (int): mass of the earth in kg
        R (float): radius of the earth in m
        rho (float): density of water (in kg/m3)
        k (float): vector of Load Love numbers until same degree as given cnm and snm.
        spacing (int, optional): Spacing of the grid in degrees. Defaults to 1.
        area (list, optional): Area for which the EWH shall be calculated. Defaults to None and the whole earth is used.
        
    Returns:
        dict: Dataset containing the EWH with keys "longitudes", "colatitudes" and "data".
            longitudes (list): List of longitudes.
            colatitudes (list): List of colatitudes.
            data (list): Numpy array containing the EWH.
    """
    if(area is None):
        # Get the pixels from the whole earth
        pixels = []
        for longitude in range(-180, 180, spacing):
            for latitude in range(-90, 90, spacing):
                pixels.append([longitude, latitude])
        pixels = np.array(pixels)
        area_weight = None
    else:
        # Get the pixels from the given area
        grid = fu.getGridfromPolygon(area, spacing)
        pixels = grid[:, :2]
        area_weight = grid[:, 2]

    # Create the vectors for the coordinates
    lamda = np.radians(pixels[:, 0])
    theta = np.radians(pixels[:, 1]+90)

    
    # Calculate the EWH (either fast or slow)
    result = calc_EWH(lamda, theta, np.array(dataset["C"]), np.array(dataset["S"]), M, R, rho, k)
    # result = fu.calc_EWH_fast(lamda, theta, np.array(dataset["C"]), np.array(dataset["S"]), M, R, rho, k)

    # Get the longitudes and colatitudes
    longitudes = []
    latitudes = []
    for pixel in pixels:
        longitudes.append(float(pixel[0]))
        latitudes.append(float(pixel[1]))
    longitudes = sorted(list(set(longitudes)))
    latitudes = sorted(list(set(latitudes)))
    
    result_matrix = np.zeros((len(latitudes), len(longitudes)))
    for index, entry in enumerate(result):
        result_matrix[int(index/len(longitudes))][index%len(longitudes)] = float(entry)
    
    new_dataset = {"longitudes": longitudes, "latitudes": latitudes, "data": result_matrix, "area_weights": area_weight}
    return(new_dataset)


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


def apply_gaussian_filtering(dataset, degree="auto", filter_radius=200):
    """
    This function is used to apply the gaussian filtering a dataset.

    Args:
        dataset (dict): Dataset containing the coefficients.
        degree (int): Maximum degree for which the factors are calculated. Defaults to "auto".
        filter_radius (int, optional): Filter radius in kilometers. Defaults to 200000.
    """
    for scalar in ["C", "S"]:
        if(degree == "auto"):
            degree = np.shape(dataset[scalar])[0]
        w = gaussian_filtering_factors(degree, filter_radius)
        for i in range(len(dataset[scalar])):
            for j in range(len(dataset[scalar][i])):
                dataset[scalar][i][j] *= w[i]
    return(dataset)
    

# Beginning of the Main Programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    
    # Importing all datasets
    # - - - - - - - - - - - -

    for index, dataset in enumerate(main.datasets):
        print(f'[Info][{index+1}/{len(main.datasets)}] Importing dataset: {dataset["name"]}')
        main.datasets[index]["data"] = import_data(dataset)
    print(f'[Info] Importing datasets finished')

    print(f'[Done] Task A')


    # Create dataset of the augmented gravity field model of each month
    # -----------------------------------------------------------------

    # Create new dataset for the output
    new_dataset = {"name": "grace_augmented", "data": []}

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
        new_dataset["data"].append({"date": itsg_grace["date"], "data": {"C": itsg_grace_matrix_c,
                                                                         "S": itsg_grace_matrix_s,
                                                                         "sigma_C": itsg_grace_matrix_sigma_c,
                                                                         "sigma_S": itsg_grace_matrix_sigma_s}})

    # Append all the datasets to the main datasets
    main.datasets.append(new_dataset)
    print(f'[Info][Done] Creating dataset of the augmented gravity field models')

    # Remove the temporary variables
    del new_dataset
    del itsg_grace_2018_matrix_c, itsg_grace_2018_matrix_s, itsg_grace_2018_matrix_sigma_c, itsg_grace_2018_matrix_sigma_s
    del itsg_grace_matrix_c, itsg_grace_matrix_s, itsg_grace_matrix_sigma_c, itsg_grace_matrix_sigma_s
    del deg1_matrix_c, deg1_matrix_s, deg1_matrix_sigma_c, deg1_matrix_sigma_s

    print(f'[Done] Task B')


    # Calculate the unfiltered equivalent water height
    # ------------------------------------------------
    
    print(f'[Info] Calculating the unfiltered equivalent water height', end="\r")
    
    selected_date = "2008-04"  # April 2008
    selected_grace = main.select_dataset(main.select_dataset(main.datasets, "name", "grace_augmented")["data"], "date", selected_date)["data"]

    # Loading the love numbers
    love_numbers = main.select_dataset(main.datasets, "name", "loadLoveNumbers_Gegout97.txt")["data"]
    
    new_dataset = apply_ewh(selected_grace, mass, radius, rho_water, love_numbers, spacing=grid_spacing, area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
    main.datasets.append({"name": f'ewh_{selected_date}',
                          "data": new_dataset["data"],
                          "axis": [new_dataset["longitudes"], new_dataset["latitudes"]]})
    del new_dataset  # Remove the temporary variable
    print(f'[Info][Done] Calculating the unfiltered equivalent water height')

    print(f'[Done] Task C')


    # Filtering of the spherical harmonic coefficients
    # ------------------------------------------------
    
    filter_radius = 300  # km

    # Create new dataset for the filtered spherical harmonic coefficients
    print(f'[Info] Creating new dataset for the filtered spherical harmonic coefficients', end="\r")
    grace_single_filtered = {"name": "gaussian_filter_coefficients"}

    # Filter the selected grace_augmented dataset for the month
    grace_single_filtered = apply_gaussian_filtering(selected_grace, filter_radius=filter_radius)
    
    main.datasets.append({"name": f'grace_{selected_date}_filtered_{filter_radius}km',
                          "data": grace_single_filtered})
    
    print(f'[Info][Done] Creating new dataset for the filtered spherical harmonic coefficients')
    
    # Create dataset with the ewh for the region of interest (filtered)
    print(f'[Info] Creating dataset with the ewh for the region of interest (filtered)', end="\r")
    new_dataset = apply_ewh(grace_single_filtered, mass, radius, rho_water, love_numbers, spacing=grid_spacing, area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
    main.datasets.append({"name": f'ewh_{selected_date}_filtered_{filter_radius}km',
                          "data": new_dataset["data"],
                          "axis": [new_dataset["longitudes"], new_dataset["latitudes"]]})
    main.datasets.append({"name": f'area_weights_{selected_date}_filtered_{filter_radius}km',
                          "data": new_dataset["area_weights"],
                          "axis": [new_dataset["longitudes"], new_dataset["latitudes"]]})
    del new_dataset  # Remove the temporary variable
    print(f'[Info][Done] Creating dataset with the ewh for the region of interest (filtered)')

    # Again, but with different filter radii
    print(f'[Info] Creating dataset with the ewh for the region of interest (filtered) for different filter radii', end="\r")
    filter_radii = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    for index, filter_radius in enumerate(filter_radii):
        print(f'[Info][{index}/{len(filter_radii)}] Creating dataset with the ewh for the region of interest (filtered) for different filter radii', end="\r")
        grace_single_filtered = apply_gaussian_filtering(selected_grace, filter_radius=filter_radius)
        new_dataset = apply_ewh(grace_single_filtered, mass, radius, rho_water, love_numbers, spacing=grid_spacing, area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
        main.datasets.append({"name": f'ewh_{selected_date}_filtered_{filter_radius}km',
                              "data": new_dataset["data"],
                              "axis": [new_dataset["longitudes"], new_dataset["latitudes"]]})
        main.datasets.append({"name": f'area_weights_{selected_date}_filtered_{filter_radius}km',
                              "data": new_dataset["area_weights"],
                              "axis": [new_dataset["longitudes"], new_dataset["latitudes"]]})
        del new_dataset  # Remove the temporary variable
    print(f'[Info][Done] Creating dataset with the ewh for the region of interest (filtered) for different filter radii')

    # Computing monthly solutions with a selected filter radius

    filter_radius = 300  # km
    
    # Create new dataset for the monthly spherical harmonic coefficients
    new_dataset = {"name": f"grace_augmented_filtered_{filter_radius}km", "data": []}

    length = len(main.select_dataset(main.datasets, "name", "grace_augmented")["data"])  # Just for the progress bar
    for index, dataset in enumerate(main.select_dataset(main.datasets, "name", "grace_augmented")["data"]):
        print(f'[Info][{index}/{length}] Computing monthly solution "{dataset["date"]}" with the selected filter radius {filter_radius} km', end="\r")  # Progress bar
        new_dataset["data"].append(apply_gaussian_filtering(dataset["data"], filter_radius=filter_radius))
    
    print(f'[Info][Done] Computing monthly solution with the selected filter radius {filter_radius} km')
        
        
    # Add the new dataset to the main datasets
    main.datasets.append(new_dataset)
    
    # Delete the temporary variables
    del new_dataset, length, index, dataset
    del filter_radii, grace_single_filtered, selected_grace

    print(f'[Done] Task D')


    # Computation of a time series of region averages
    # -----------------------------------------------

    # Create new dataset for the monthly solutions of equivalent water height and area weights
    new_dataset_ehw = {"name": f"monthly_ewh_filtered_{filter_radius}km", "data": []}
    new_dataset_area_weights = {"name": f"monthly_area_weights_filtered_{filter_radius}km", "data": []}

    length = len(main.select_dataset(main.datasets, "name", "grace_augmented_filtered_{filter_radius}km")["data"])  # Just for the progress bar
    for index, dataset in enumerate(main.select_dataset(main.datasets, "name", "grace_augmented_filtered_{filter_radius}km")["data"]):
        print(f'[Info][{index}/{length}] Computing monthly solution "{dataset["date"]}" with the selected filter radius {filter_radius} km', end="\r")  # Progress bar
        ewh = apply_ewh(dataset, mass, radius, rho_water, love_numbers, spacing=grid_spacing, area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
        new_dataset_ehw["data"].append({"date": dataset["date"], "data": ewh["data"], "axis": [ewh["longitudes"], ewh["latitudes"]]})
        new_dataset_area_weights["data"].append({"date": dataset["date"], "data": ewh["area_weights"], "axis": [ewh["longitudes"], ewh["latitudes"]]})

    # Add the new datasets to the main datasets
    main.datasets.append(new_dataset_ehw)
    main.datasets.append(new_dataset_area_weights)

    # Delete the temporary variables
    del new_dataset_ehw, new_dataset_area_weights, length, index, dataset, ewh










    # Interpolation of missing GRACE months
    # -----------------------------------------------

    # grace_missing_months = fu.interp_missing_months()

    # Estimation of the linear mass trend
    # -----------------------------------------------



    # Attempt to save all datasets
    
    for dataset in main.datasets:
        export_data(dataset)
