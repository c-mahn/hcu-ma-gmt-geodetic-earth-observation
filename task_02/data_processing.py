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

gravity_constant = 3.986005000e+14   # m^3 / (kg * s^2)
normal_gravity = 9.80665           # m / s^2
radius = 6.378137000e+06   # m
mass = 5.9722*10**24     # kg
degree_n = 80
order_m = degree_n/2
rho_grad = 180/np.pi
rho_water = 1000              # kg / m^3
grid_spacing = 1                 # degree


# Other Variables
# -----------------------------------------------------------------------------

first_year = 2003
last_year = 2016
months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
filter_radii = [100, 200, 300, 400, 500]
filter_radius = 300  # km
selected_date = "2008-04"  # April 2008

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
    ignore_until = "end_of_head"
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
            if(line[0] == "gfc" and int(line[1]) <= degree_n):
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
        data = np.array(data)
    else:
        print(f'[Error] The dataset {dataset["name"]} has an unknown type')
        return(None)
    return(data)


def export_data(dataset):
    try:
        if(type(dataset["data"]) is list):
            for subdataset in dataset["data"]:
                with open(os.path.join("output", f'{subdataset["name"]}.csv'), 'w') as f:
                    for index, value in enumerate(subdataset["ewh"]):
                        f.write(f'{subdataset["lamda"][index]};{subdataset["theta"][index]};{value}\n')
    except:
        pass
    try:
        with open(os.path.join("output", f'{dataset["name"]}.csv'), 'w') as f:
            for index, value in enumerate(dataset["ewh"]):
                f.write(f'{dataset["lamda"][index]};{dataset["theta"][index]};{value}\n')
    except:
        pass


def assemble_matrix(data, value_index, coord_indices=["L", "M"]):
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
    matrix = np.zeros((dimension_x+1, dimension_y+1))
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
        # Calculate the legendre functions
        P = fu.legendreFunctions(theta_i, max_degree)
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
        area_weight = 1
    else:
        # Get the pixels from the given area
        grid = fu.getGridfromPolygon(area, spacing)
        pixels = grid[:, :2]
        area_weight = grid[:, 2]

    # Create the vectors for the coordinates
    lamda = np.radians(pixels[:, 0])
    theta = np.radians(90-pixels[:, 1])

    # Cut the k vector to the same length as the cnm and snm
    # Needed for the fast version (in the slow version the other values are ignored)
    k = k[0:np.shape(dataset["C"])[0]]

    # Calculate the EWH (either fast or slow)
    # result = calc_EWH(lamda, theta, np.array(dataset["C"]), np.array(dataset["S"]), M, R, rho, k)
    result = fu.calc_EWH_fast(lamda, theta, dataset["C"], dataset["S"], M, R, rho, k)

    # Create the new dataset and return it
    new_dataset = {"lamda": pixels[:, 0], "theta": pixels[:,
                                                          1], "ewh": result, "area_weights": area_weight}
    return(new_dataset)


def interpolate_missing_data(values):
    """
    This function is used to interpolate the missing data in the datasets.
    """
    # create date vector
    startyr = 2003
    endyear = 2016
    mn = np.array(range(1, 13))
    yr = np.array(range(startyr, endyear+1))

    # convert to decimal years
    dec_yr = []
    for year in yr:
        for month in mn:
            dec_yr.append(year+((month-1)/12))

    dec_yr = np.array(dec_yr)

    # interpolate missing data
    for i in range(np.shape(values)[0]):
        if(values[i] == 0):
            values[i] = np.mean([float(values[i-1]), float(values[i+1])])
    return(dec_yr, values)


def gaussian_filtering_factors(degree, filter_radius=500):
    """
    This function is used to calculate the gaussian filtering factors.

    Args:
        degree (int): Maximum degree for which the factors are calculated.
        filter_radius (int, optional): Filter radius in kilometers. Defaults to 500.
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


def apply_gaussian_filtering(dataset, filter_radius=500):
    """
    This function is used to apply the gaussian filtering to a dataset.

    Args:
        dataset (dict): Dataset containing the coefficients.
        degree (int): Maximum degree for which the factors are calculated. Defaults to "auto".
        filter_radius (int, optional): Filter radius in kilometers. Defaults to 200000.
    """
    for scalar in ["C", "S"]:
        dataset[scalar] = dataset[scalar].tolist()
        degree = np.shape(dataset[scalar])[0]
        w = gaussian_filtering_factors(degree, filter_radius)
        for n in range(len(dataset[scalar])):
            for m in range(len(dataset[scalar][n])):
                dataset[scalar][n][m] *= w[n]
        dataset[scalar] = np.array(dataset[scalar])
    return(dataset)


def plot_xy(diagrams, names=["measurement"], title="automatic", plot="show", axis={"x": "x", "y": "y"}, log={"x": False, "y": False}):
    """
    This function is used to plot graphs in the x-y-plane.

    Args:
        graphs (list): List of graphs to plot. Defaults to [{"x": [0, 1], "y": [0, 1]}] as an example.
        names (list, optional): List of names for the graphs. Defaults to ["measurement"].
        title (str, optional): Title of the plot. Defaults to "automatic". When set to "automatic" the title will be the same as the first name.
        plot (str, optional): Defines if the plot should be shown or saved. Defaults to "show".
    """
    if(log["x"]):
        plt.xscale("log")
    if(log["y"]):
        plt.yscale("log")
    if(title == "automatic"):
        title = names[0]
    for item in diagrams:
        plt.plot(item["x"], item["y"])
    plt.legend(names)
    plt.grid()
    plt.xlabel(axis["x"])
    plt.ylabel(axis["y"])
    plt.title(title)
    if(plot == "show"):
        plt.show()
    elif(plot == "save"):
        plt.savefig(title + ".png")
    plt.cla()


# Beginning of the Main Programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Importing all datasets
    # - - - - - - - - - - - -

    for index, dataset in enumerate(main.datasets):
        print(
            f'[Info][{index+1}/{len(main.datasets)}] Importing dataset: {dataset["name"]}')
        main.datasets[index]["data"] = import_data(dataset)
    print(f'[Info] Importing datasets finished')

    print(f'[Done] Task A')

    # Create dataset of the augmented gravity field model of each month
    # -----------------------------------------------------------------

    # Create new dataset for the output
    new_dataset = {"name": "grace_augmented", "data": []}

    # Load the static grace observation model and assemble it into matrices
    itsg_grace_static = main.select_dataset(
        main.datasets, "name", "ITSG-Grace2018s.gfc")
    for variable in ["C", "S", "sigma_C", "sigma_S"]:
        itsg_grace_static[variable] = assemble_matrix(
            itsg_grace_static["data"], value_index=variable)

    length = len(main.select_dataset(main.datasets, "name",
                 "ITSG-Grace")["data"])  # Just for the progress bar
    for index, itsg_grace in enumerate(main.select_dataset(main.datasets, "name", "ITSG-Grace")["data"]):
        # Progress bar
        print(
            f'[Info][{index+1}/{length}] Creating augmented grace-dataset ({itsg_grace["date"]})', end="\r")
        new_itsg_grace_dataset = itsg_grace

        # Load the corresponding deg1 coefficients
        deg1_dataset = []
        for dataset in main.select_dataset(main.datasets, "name", "deg1")["data"]:
            if(dataset["date"] == itsg_grace["date"]):  # Selecting the correct dataset
                deg1_dataset = dataset
                break

        for variable in ["C", "S", "sigma_C", "sigma_S"]:
            # Assemble monthly grace observations into matrices
            new_itsg_grace_dataset[variable] = assemble_matrix(new_itsg_grace_dataset["data"], value_index=variable)

            # - static grace model
            new_itsg_grace_dataset[variable] = matrix_math(new_itsg_grace_dataset[variable], itsg_grace_static[variable], operator="-")

            # Assemble the deg1 coefficients into matrices
            deg1_dataset[variable] = assemble_matrix(deg1_dataset["data"], value_index=variable)

            # + monthly grace coefficients
            new_itsg_grace_dataset[variable] = matrix_math(new_itsg_grace_dataset[variable], deg1_dataset[variable], operator="+")

        # Append the monthly augmented grace dataset
        new_dataset["data"].append({"date": itsg_grace["date"], "data": {"C": new_itsg_grace_dataset["C"],
                                                                         "S": new_itsg_grace_dataset["S"],
                                                                         "sigma_C": new_itsg_grace_dataset["sigma_C"],
                                                                         "sigma_S": new_itsg_grace_dataset["sigma_S"]}})

    # Append all the datasets to the main datasets
    main.datasets.append(new_dataset)
    print(f'[Info][Done] Creating dataset of the augmented gravity field models')

    # Remove the temporary variables
    del new_dataset, itsg_grace_static, deg1_dataset, new_itsg_grace_dataset

    print(f'[Done] Task B')

    # Calculate the unfiltered equivalent water height
    # ------------------------------------------------

    print(f'[Info] Calculating the unfiltered equivalent water height', end="\r")
    selected_grace = main.select_dataset(main.select_dataset(main.datasets, "name", "grace_augmented")["data"], "date", selected_date)["data"]

    # Loading the love numbers
    love_numbers = main.select_dataset(main.datasets, "name", "loadLoveNumbers_Gegout97.txt")["data"]

    new_dataset = apply_ewh(selected_grace, mass, radius, rho_water, love_numbers, spacing=grid_spacing,
                            area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
    new_dataset["name"] = f'ewh_{selected_date}'
    main.datasets.append(new_dataset)
    del new_dataset  # Remove the temporary variable
    print(f'[Info][Done] Calculating the unfiltered equivalent water height')

    print(f'[Done] Task C')

    # Filtering of the spherical harmonic coefficients
    # ------------------------------------------------

    # Create new dataset for the filtered spherical harmonic coefficients
    grace_single_filtered = apply_gaussian_filtering(selected_grace, filter_radius=filter_radius)
    main.datasets.append({"name": f'grace_{selected_date}_filtered_{filter_radius}km',
                          "data": grace_single_filtered})

    print(f'[Info][Done] Creating new dataset for the filtered spherical harmonic coefficients')

    # Create dataset with the ewh for the region of interest (filtered)
    print(f'[Info] Creating dataset with the ewh for the region of interest (filtered)', end="\r")
    new_dataset = apply_ewh(grace_single_filtered, mass, radius, rho_water, love_numbers, spacing=grid_spacing,
                            area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
    new_dataset["name"] = f'ewh_{selected_date}_filtered_{filter_radius}km'
    main.datasets.append(new_dataset)
    del new_dataset, grace_single_filtered  # Remove the temporary variable
    print(f'[Info][Done] Creating dataset with the ewh for the region of interest (filtered)')

    # Again, but with different filter radii
    temp_list_for_plotting = []
    print(f'[Info] Creating datasets with the ewh for the region of interest (filtered) for different filter radii', end="\r")
    for index, radius_i in enumerate(filter_radii):
        print(f'[Info][{index}/{len(filter_radii)}] Creating dataset with the ewh for the region of interest (filtered) for different filter radii', end="\r")
        grace_single_filtered = apply_gaussian_filtering(selected_grace, filter_radius=radius_i)
        new_dataset = apply_ewh(grace_single_filtered, mass, radius, rho_water, love_numbers, spacing=grid_spacing,
                                area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"])
        new_dataset["name"] = f'ewh_{selected_date}_filtered_for_{radius_i}km'
        new_dataset_2 = grace_single_filtered
        new_dataset_2["name"] = f'grace_{selected_date}_filtered_for_{radius_i}km'
        main.datasets.append(new_dataset)
        main.datasets.append(new_dataset_2)
        temp_list_for_plotting.append(new_dataset_2)
    del new_dataset, new_dataset_2, radius_i  # Remove the temporary variable
    print(f'[Info][Done] Creating dataset with the ewh for the region of interest (filtered) for different filter radii')

    # Plot of the gaussian filtering factors

    # List of all the gaussian filtering factors for different filter radii
    gaussian_factor_sets = []
    for radius_i in filter_radii:
        gaussian_factor_sets.append(gaussian_filtering_factors(degree_n, radius_i))

    graphs = []  # List of all the graphs, that will be plotted in the diagram
    for factors in gaussian_factor_sets:
        x_values = []
        y_values = []
        for index, factor in enumerate(factors):
            x_values.append(index)
            y_values.append(float(factor))
        graphs.append({"x": x_values, "y": y_values})

    plot_xy(graphs,
            names=filter_radii,
            title=f'Gaussian filtering factors for different filter radii',
            axis={"x": "Degree", "y": "Factor"},
            plot="save")

    # Remove the temporary variable
    del gaussian_factor_sets, index, radius_i, graphs, x_values, y_values, factors, factor

    # Plot of the degree variances for different filter radii

    graphs = []  # List of all the graphs, that will be plotted in the diagram
    names = []  # List of all the names of the graphs

    degree_variances = []  # List of all the degree variances for different filter radii
    for index, radius_i in enumerate(filter_radii):
        dataset = temp_list_for_plotting[index]
        graph_c = {"x": [], "y": []}
        graph_s = {"x": [], "y": []}
        c = dataset["C"]
        s = dataset["S"]
        for n, line in enumerate(c):
            variance_sum_c = 0
            variance_sum_s = 0
            for m, entry in enumerate(c[n]):
                variance_sum_c += c[n][m]**2
                variance_sum_s += s[n][m]**2
            graph_c["x"].append(n)
            graph_c["y"].append(variance_sum_c)
            graph_s["x"].append(n)
            graph_s["y"].append(variance_sum_s)
        graphs.append(graph_c)
        graphs.append(graph_s)
        names.append(f'C_{radius_i}km')
        names.append(f'S_{radius_i}km')

    for index, graph in enumerate(graphs):
        plot_xy([graph],
                names=f'{names[index]}',
                title=f'Degree variances for different filter radii ({names[index]})',
                axis={"x": "Degree", "y": "Variance"},
                log={"x": False, "y": True},
                plot="save")

    # Computing monthly solutions with a selected filter radius
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Create new dataset for the monthly spherical harmonic coefficients
    new_dataset = {"name": f'monthly_grace_coefficients_filtered_{filter_radius}km', "data": []}

    length = len(main.select_dataset(main.datasets, "name", "grace_augmented")["data"])  # Just for the progress bar
    for index, dataset in enumerate(main.select_dataset(main.datasets, "name", "grace_augmented")["data"]):
        print(f'[Info][{index}/{length}] Computing monthly datasets "{dataset["date"]}" with filter radius of {filter_radius} km', end="\r")
        new_dataset["data"].append(apply_gaussian_filtering(
            dataset["data"], filter_radius=filter_radius))
        # Add the date to the dataset
        new_dataset["data"][-1]["date"] = dataset["date"]
    main.datasets.append(new_dataset)

    print(f'[Info][Done] Computing monthly datasets with the filtered spherical harmonic coefficients with the selected filter radius {filter_radius} km')

    # Delete the temporary variables
    del new_dataset, length, index, dataset
    del grace_single_filtered, selected_grace

    print(f'[Done] Task D')

    # Computation of a time series of region averages
    # -----------------------------------------------

    # Create new dataset for the monthly equivalent water height

    new_dataset = {"name": f'monthly_ewh_filtered_{filter_radius}km_region_bounding_box', "data": []}
    new_dataset_2 = {"name": f'monthly_ewh_filtered_{filter_radius}km_region_polygon', "data": []}

    length = len(main.select_dataset(main.datasets, "name", f'monthly_grace_coefficients_filtered_{filter_radius}km')["data"])  # Just for the progress bar
    for index, dataset in enumerate(main.select_dataset(main.datasets, "name", f'monthly_grace_coefficients_filtered_{filter_radius}km')["data"]):
        # Progress bar
        print(f'[Info][{index}/{length}] Computing monthly solution "{dataset["date"]}" with filter radius of {filter_radius} km', end="\r")
        new_dataset["data"].append(apply_ewh(dataset, mass, radius, rho_water, love_numbers, spacing=grid_spacing,
                                             area=main.select_dataset(main.datasets, "name", "region_bounding_box.txt")["data"]))
        new_dataset["data"][-1]["name"] = f'ewh_{dataset["date"]}_filtered_{filter_radius}km_region_bounding_box'
        # Add the date to the dataset
        new_dataset["data"][-1]["date"] = dataset["date"]
        new_dataset_2["data"].append(apply_ewh(dataset, mass, radius, rho_water, love_numbers, spacing=grid_spacing,
                                               area=main.select_dataset(main.datasets, "name", "region_polygon.txt")["data"]))
        new_dataset_2["data"][-1]["name"] = f'ewh_{dataset["date"]}_filtered_{filter_radius}km_region_polygon'
        # Add the date to the dataset
        new_dataset_2["data"][-1]["date"] = dataset["date"]

        # Apply the area weighting
        new_dataset["data"][-1]["ewh"] *= new_dataset["data"][-1]["area_weights"]*radius**2
        new_dataset_2["data"][-1]["ewh"] *= new_dataset_2["data"][-1]["area_weights"]*radius**2
        new_dataset["data"][-1]["ewh"] /= 1e9
        new_dataset_2["data"][-1]["ewh"] /= 1e9

        # Compute mean
        new_dataset["data"][-1]["mean"] = np.mean(new_dataset["data"][-1]["ewh"])
        new_dataset_2["data"][-1]["mean"] = np.mean(new_dataset_2["data"][-1]["ewh"])

    main.datasets.append(new_dataset)
    main.datasets.append(new_dataset_2)

    print(f'[Info][Done] Computing monthly solutions of equivalent water height with the selected filter radius {filter_radius} km')

    print(f'Area of region bounding box: {np.sum(new_dataset["data"][-1]["area_weights"]*radius**2)/1000000:.2f} km??')
    print(f'Area of region polygon: {np.sum(new_dataset_2["data"][-1]["area_weights"]*radius**2)/1000000:.2f} km??')

    # Delete the temporary variables
    del length, index, dataset, new_dataset, new_dataset_2

    print(f'[Done] Task E')

    # Interpolation of missing GRACE months
    # -----------------------------------------------

    # Create new dataset for the collected monthly equivalent water height means
    new_dataset = {"name": f'collection_of_monthly_ewh_means_f{filter_radius}', "data": []}
    for dataset in main.select_dataset(main.datasets, "name", f'monthly_ewh_filtered_{filter_radius}km_region_polygon')["data"]:
        new_dataset["data"].append({"date": dataset["date"], "mean": dataset["mean"]})
    main.datasets.append(new_dataset)

    # Delete the temporary variables
    del new_dataset

    # Interpolate the missing months
    temp_vector = np.zeros(len(months)*(last_year-first_year+1))  # *np.nan
    graph_means = {"x": [], "y": []}
    for dataset in main.select_dataset(main.datasets, "name", f'collection_of_monthly_ewh_means_f{filter_radius}')["data"]:
        month = int(dataset["date"].split("-")[1])
        year = int(dataset["date"].split("-")[0])
        temp_vector[(year-first_year)*len(months)+month-1] = dataset["mean"]
    dates, means = interpolate_missing_data(temp_vector)

    # Save the interpolated monthly means
    new_dataset = {"name": f'interpolated_monthly_ewh_means_f{filter_radius}', "ewh": means, "dates": dates}

    # Plot the uninterpolated monthly means
    graph_means["y"] = temp_vector
    graph_means["x"] = dates
    plot_xy([graph_means],
            names=["Monthly means"],
            title=f'Monthly means of the equivalent water height with a filter radius of {filter_radius} km',
            axis={"x": "Time [years]", "y": "Equivalent water height [m]"},
            plot="save")

    main.datasets.append(new_dataset)

    # Plot the interpolated monthly means
    graph_means["y"] = means
    graph_means["x"] = dates
    plot_xy([graph_means],
            names=["Monthly means"],
            title=f'Interpolated monthly means of the equivalent water height with a filter radius of {filter_radius} km',
            axis={"x": "Time [years]", "y": "Equivalent water height [m]"},
            plot="save")

    # Delete the temporary variables
    del temp_vector, dates, means, graph_means, new_dataset

    print(f'[Done] Task F')

    # Estimation of the linear mass trend
    # -----------------------------------------------

    # Calculate the linear trend
    timeseries = main.select_dataset(main.datasets, "name", f'interpolated_monthly_ewh_means_f{filter_radius}')["ewh"]
    time = main.select_dataset(main.datasets, "name", f'interpolated_monthly_ewh_means_f{filter_radius}')["dates"]
    timeseries = np.array(timeseries)
    time = np.array(time)
    slope = (len(timeseries) * np.sum(time*timeseries) - np.sum(time) * np.sum(timeseries)) / (len(time)*np.sum(time*time) - np.sum(time) ** 2)

    print(f'[Info] The linear trend of the equivalent water height is {slope:.2f} m/yr')

    # Convert to gigatonnes per year
    gigatonnes = slope * mass * rho_water / 1e9**3

    print(f'[Info] The linear trend of the equivalent water height is {gigatonnes:.2f} Gt/yr')
    
    # Plot the linear trend
    graphs = [{"x": time, "y": (timeseries-timeseries[0]) * mass * rho_water / 1e9**3}]
    graphs.append({"x": time, "y": (time-first_year)*slope * mass * rho_water / 1e9**3})
    plot_xy(graphs,
            names=["Monthly means", "Linear trend"],
            title=f'Linear trend of the equivalent water height with a filter radius of {filter_radius} km',
            axis={"x": "Time [years]", "y": "Equivalent water height [Gt]"},
            plot="save")

    # Delete the temporary variables
    del timeseries, time, slope, gigatonnes

    print(f'[Done] Task G')

    # Post-processing
    # -----------------------------------------------

    # Attempt to save all datasets
    for index, dataset in enumerate(main.datasets):
        print(f'[Info][{index+1}/{len(main.datasets)}] Exporting dataset "{dataset["name"]}"')
        export_data(dataset)
