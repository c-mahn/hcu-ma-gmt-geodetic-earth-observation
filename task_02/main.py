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

# These are the datasets that will be used for the computation.

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
