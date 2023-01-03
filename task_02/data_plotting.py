# Data Plotting
# #############################################################################

# This python script plots the spherical harmonics using the commandline-tool
# gmt.

# Authors:
# Christopher Mahn
# Silas Teske
# Joshua Wolf

# #############################################################################

# Import of Libraries
# -----------------------------------------------------------------------------

import os
import shutil


# Functions
# -----------------------------------------------------------------------------

def plot_sherical_harmonics(plot_name="plot",
                            img_type="png",
                            title="Title",
                            subtitle="Subtitle",
                            nc_filename="input.nc",
                            colorbar_setting='-Dx0c/-1.5c+w12c/0.25c+h -B20+l"[m]" -V',
                            show_plot=False):
    """
    This function plots spherical harmonics using the commandline-tool gmt.
    GMT is a free and open-source command-line tool for plotting and processing
    data. It is available for Windows, Linux and Mac. It can be downloaded here:
    https://www.generic-mapping-tools.org/download/
    
    Args:
        plot_name (str): Name of the plot
        img_type (str): Type of the image (eg: "png")
        title (str): Title of the plot
        subtitle (str): Subtitle of the plot
        nc_filename (str): Name of the netCDF-file
        colorbar_setting (str): Settings for the colorbar
        show_plot (bool): If True, the plot will be shown after it has been saved
    """
    os.system(f'gmt begin {plot_name} {img_type}')  # Begin the plot
    os.system(f'gmt basemap -JN12c -R-180/180/-90/90 -B+t"{title}" -B')  # Set the basemap
    os.system(f'gmt grd2cpt {nc_filename} -Ccyclic.cpt -Sh -E100 -Z -V')  # Set the colorbar
    os.system(f'gmt grdimage {nc_filename} -I+a315+ne0.2')  # Plot the data
    os.system(f'gmt coast -W.1,120/120/120 -B')  # Plot the coast
    os.system(f'gmt text -F+cBL+t"{subtitle}" -N -D2.3c/-1c')  # Plot the subtitle
    os.system(f'gmt colorbar {colorbar_setting}')  # Plot the colorbar
    if(show_plot):  # Show the plot, if specified
        os.system(f'gmt end show')  # Show the plot
    else:  # Otherwise, just save the plot
        os.system(f'gmt end')  # Save the plot
    
    # Move file to the plots folder
    shutil.move(f'{plot_name}.{img_type}', f'./plots/{plot_name}.{img_type}')


# Classes
# -----------------------------------------------------------------------------


# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # Plotting of the Gravity-Anomalies and Geoid-Height
    plot_sherical_harmonics(plot_name="geoid_height",
                            title="Geoid Height",
                            subtitle="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            nc_filename="geoid_height.nc",
                            colorbar_setting='-Dx0c/-1.5c+w12c/0.25c+h -B20+l"[m]" -V')

    plot_sherical_harmonics(plot_name="grav_anom_satellite",
                            title="Gravitational Anomalies Satellite",
                            subtitle="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            nc_filename="grav_anom_satellite.nc",
                            colorbar_setting='-Dx0c/-1.5c+w12c/0.25c+h -B0.05+l"[mm]" -V')

    plot_sherical_harmonics(plot_name="grav_anom_surface",
                            title="Gravitational Anomalies Surface",
                            subtitle="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            nc_filename="grav_anom_surface.nc",
                            colorbar_setting='-Dx0c/-1.5c+w12c/0.25c+h -B0.1+l"[mm]" -V')
