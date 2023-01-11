# Main-Script
# #############################################################################

# This python script automatically launches all other python scripts in the
# right order and computes the entire task.

# Authors:
# Christopher Mahn
# Silas Teske
# Joshua Wolf

# #############################################################################

# Import of Libraries
# -----------------------------------------------------------------------------

# import math as m
# import numpy as np
# import string as st
# import random as r
# import re
import os
import platform


# Project-Settings
# -----------------------------------------------------------------------------

# These settings affect how the executed scripts below will compute the data.
# Changing these values may increase execution-time significantly or allowes to
# change the computed input or output.
folder_data = "data"

# These are the datasets that will be used for the computation.
datasets = [{"name":"ITSG-Grace", "type": "gfc", "is_folder": True},
            {"name":"deg1", "type": "gfc", "is_folder": True},
            {"name": "ITSG-Grace2018s.gfc", "type": "gfc", "is_folder": False},
            {"name": "loadLoveNumbers_Gegout97.txt", "type": "single_column", "is_folder": False},
            {"name": "region_bounding_box.txt", "type": "csv", "delimiter": ",", "is_folder": False},
            {"name": "region_polygon.txt", "type": "csv", "delimiter": ",", "is_folder": False}]


# Functions
# -----------------------------------------------------------------------------

def __run_script(script_name):
    """
    This function executes python scripts via the command line.

    Args:
        script_name (str): name of the python script (eg: "demo.py")
    """
    if(platform.system() == "Linux"):
        print(f'[INFO] Executing "{script_name}" as Linux-User')
        os.system(f'python3 {script_name}')  # Run on Linux
    elif(platform.system() == "Windows"):
        user = os.environ.get('USERNAME')
        print(f'[INFO] Executing "{script_name}" as Windows-User "{user}"')
        os.system(f'C:/Users/{user}/anaconda3/python.exe {script_name}')  # Run on Windows


def terminate():
    """
    This function terminates the program.
    """
    print("[INFO] The program has been terminated.")
    exit()


def select_dataset(datasets, key, value):
    """
    This function selects a dataset from the list of datasets.

    Args:
        datasets (list): list of datasets
        key (str): key of the dataset
        value (str): value of the dataset

    Returns:
        dict: dataset
    """
    for dataset in datasets:
        if(dataset[key] == value):
            return dataset
    return None


# Classes
# -----------------------------------------------------------------------------


# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    os.makedirs("processed_data", exist_ok=True)
    __run_script("data_processing.py")
    os.makedirs("plots", exist_ok=True)
    # __run_script("data_plotting.py")
    print("[INFO] All calculations have ended.")
