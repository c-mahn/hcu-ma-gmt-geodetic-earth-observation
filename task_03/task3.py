# Geodetic Earth Observation - Task 03
# #############################################################################

# Authors:
# Christopher Mahn
# Silas Teske
# Joshua Wolf

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
import numpy as np
# import math as m
# import sys
import os
import shutil
# from scipy.fft import fft, fftfreq
# from scipy import signal
# import functions as fu
# import json
import netCDF4 as nc
import functions as fu

# #############################################################################


# Constants
# -----------------------------------------------------------------------------

g                   = 9.807                 # Gravitational acceleration [m/s^2]
omega               = (2*np.pi)/(24*60*60)  # Rotational velocity of the Earth
radius              = 6.378137000e+06       # m
rho_grad            = np.pi/180             # kg/m^3

# Functions
# -----------------------------------------------------------------------------

def run_gmt(input_file_name="test_input",
            output_file_name=None,
            sample_and_cut=False,
            img_type="png",
            grid_resolution="2m",
            map_projection="B-130/65/45/65/18c",
            region ="-145/-110/45/65",
            color_palette="haxby",
            color_settings="-1/1/0.001",
            title="No Title",
            subtitle=None,
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'):
    """
    This function plots spherical harmonics using the commandline-tool gmt.
    GMT is a free and open-source command-line tool for plotting and processing
    geodata. It is available for Windows, Linux and Mac. It can be downloaded here:
    https://www.generic-mapping-tools.org/download/
    
    Args:
        file_name (str): The name of the file to be plotted.
        grid_resolution (str): The resolution of the grid to be plotted. In degrees.
        file_name_poly (str): The name of the file containing the polygon defining the area.
        map_projection (str): The map projection to be used. See GMT documentation for more information.
        region (str): The region to be plotted.
        color_palette (str): The color palette to be used. See GMT documentation for more information.
        title (str): The title of the plot.
        subtitle (str): The subtitle of the plot.
        editors (str): The names of the editors of the plot.
        colorbar_setting (str): The settings for the colorbar. See GMT documentation for more information.
        show_plot (bool): If True, the plot is shown. If False, the plot is only saved to a png file.
    """
    if output_file_name is None:
        output_file_name = input_file_name
    
    if subtitle is None:
        subtitle = output_file_name

    command = ""
    command += f'gmt gmtset FORMAT_GEO_MAP ddd'                         # Set the format of the map
    command += f' && gmt begin {output_file_name} {img_type}'           # Start the plot
    if(sample_and_cut is True):
        command += f' && gmt grdsample ./data/{input_file_name}.nc -G./data/{input_file_name}_{grid_resolution}.nc -I{grid_resolution}'   # Sample the grid
        command += f' && gmt grdcut ./data/{input_file_name}_{grid_resolution}.nc -G./data/{output_file_name}.nc -R{region}'                  # Cut the grid
    command += f' && gmt grdinfo ./data/{output_file_name}.nc'                                                         # Get the information of the grid
    command += f' && gmt makecpt -C{color_palette} -T{color_settings} -Z'                                                   # Create the color palette
    command += f' && gmt grd2cpt ./data/{output_file_name}.nc -C{color_palette} -Z'                                    # Apply the color palette
    command += f' && gmt grdimage -J{map_projection} -R{region} ./data/{output_file_name}.nc -Q'                       # Plot the grid
    command += f' && gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"{title}" -W0.25p,80/80/80 -Df -N1/1.25p,black -V '              # Plot the coastlines and the title
    command += f' && gmt text -F+cBL+t"{subtitle}" -N -D6.65c/-1c'      # Plot the subtitle
    command += f' && gmt text -F+cBL+t"{editors}" -N -D5.15c/-1.5c'     # Plot the editors
    command += f' && gmt colorbar {colorbar_settings}'                  # Plot the colorbar
    command += f' && gmt end'                                           # Save the plot                 

    os.system(command)
    
    # Move file to the plots folder
    shutil.move(os.path.join(f'{output_file_name}.{img_type}'), os.path.join('plots', f'{output_file_name}.{img_type}'))

def calc_velocity(mdt, lon, lat, grid_spacing=0.5):
    u = np.zeros((len(lat), len(lon)))
    v = np.zeros((len(lat), len(lon)))
    for index_phi, phi in enumerate(lat):
        f = 2*omega*np.sin(phi)
        for index_lam, lam in enumerate(lon):
            if index_lam == 0:
                continue
            elif index_lam == len(lon)-1:
                continue
            elif index_phi == 0:
                continue
            elif index_phi == len(lat)-1:
                continue
            else:
                dx = radius * np.cos(phi)*(grid_spacing*rho_grad)
                dy = radius * (grid_spacing*rho_grad)
                dhx = mdt[index_phi][index_lam + 1] - mdt[index_phi][index_lam - 1]
                dhy = mdt[index_phi + 1][index_lam] - mdt[index_phi - 1][index_lam]
                u[index_phi][index_lam] = -(g/f)*(dhy/(2*dy))
                v[index_phi][index_lam] = (g/f)*(dhx/(2*dx))
    return u, v



# Classes
# -----------------------------------------------------------------------------


# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Resampling and cutting the grid-data

    run_gmt(input_file_name="DTU15MDT_2min.mdt",
            output_file_name="DTU15MDT_30min_cut.mdt",
            sample_and_cut=True,
            img_type="png",
            grid_resolution="30m",
            map_projection="B-130/65/45/65/18c",
            region ="-145/-110/45/65",
            color_palette="haxby",
            color_settings="-1/1/0.001",
            title="British Columbia, Canada",
            subtitle="DTU15MDT_30min_cut.mdt.nc",
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"MDT [m]" -V' 
            )

    # Importing the grid-data with netCDF4 and saving it with the function save_grid

    data_mdt = nc.Dataset(os.path.join(f'./data/DTU15MDT_30min_cut.mdt.nc'))
    lon = data_mdt.variables["lon"][:]
    lat = data_mdt.variables["lat"][:]
    mdt = data_mdt.variables["z"][:]
    fu.save_grid(os.path.join(f'./data/DTU15MDT_30min_cut_pyout.mdt.nc'), mdt, lon, lat)

    # Plotting the grid-data saved with the function save_grid

    run_gmt(input_file_name="DTU15MDT_30min_cut_pyout.mdt",
            grid_resolution="30m",
            map_projection="B-130/65/45/65/18c",
            region ="-145/-110/45/65",
            color_settings="-1/1/0.001",
            title="British Columbia, Canada",
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"MDT [m]" -V'
            )
    
    # Calculation of geostrophic currents

    u, v = calc_velocity(mdt, lon, lat, grid_spacing=0.5)

    # Saving the geostrophic currents

    fu.save_grid(os.path.join(f'./data/velocity_components_east.mdt.nc'), u, lon, lat)
    fu.save_grid(os.path.join(f'./data/velocity_components_north.mdt.nc'), v, lon, lat)

    # Plotting the geostrophic currents

    run_gmt(input_file_name="velocity_components_east.mdt",
            grid_resolution="30m",
            map_projection="B-130/65/45/65/18c",
            region ="-145/-110/45/65",
            color_settings="-20/20/5",
            title="British Columbia, Canada",
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.2+l"Velocity [m/s]" -V'
            )
    
    run_gmt(input_file_name="velocity_components_north.mdt",
            grid_resolution="30m",
            map_projection="B-130/65/45/65/18c",
            region ="-145/-110/45/65",
            color_settings="-25/25/2.5",
            title="British Columbia, Canada",
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.2+l"Velocity [m/s]" -V'
            )