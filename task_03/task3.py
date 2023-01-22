# Geodetic Earth Observation - Task 03
# #############################################################################

# Authors:
# Silas Teske
# Joshua Wolf
# Christopher Mahn

# #############################################################################

# Import of Libraries
# -----------------------------------------------------------------------------

# import main
# import string as st
# import random as r
# import re
# from turtle import position
# import matplotlib.pyplot as plt
# from scipy import interpolate
# import numpy as np
# import math as m
# import sys
import os
# from scipy.fft import fft, fftfreq
# from scipy import signal
# import functions as fu
# import json
import pymgt

new_grid = pygmt.grdsample(grid='file=DTU15MDT_2min.mdt.nc', resolution=30m, region=[-145, -110, 45, 65], newgridname='earth_relief_30min.mdt.nc')