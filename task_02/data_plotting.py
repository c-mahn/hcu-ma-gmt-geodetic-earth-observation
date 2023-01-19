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
import main


# Functions
# -----------------------------------------------------------------------------

def plot_sherical_harmonics(file_name="test",
                            img_type="png",
                            grid_resolution="0.5",
                            file_name_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="B-130/65/45/65/18c",
                            region ="-145/-110/45/65",
                            color_palette="haxby",
                            title="No Title",
                            subtitle="No Subtitle",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V',
                            show_plot=False):
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

    """    os.system(f'gmt gmtset FORMAT_GEO_MAP ddd')  # Set the format of the map
    os.system(f'set ascii_file={file_name}.txt')  # Set the ascii file
    os.system(f'set grid_file={file_name}.grd')  # Set the grid file
    os.system(f'gmt begin {file_name} {img_type}')  # Start the plot
    os.system(f'gmt xyz2grd %ascii_file% -R{region} -r -I{grid_resolution} -G%grid_file% -V')  # Convert the ascii file to a grid file
    os.system(f'gmt grd2cpt %grid_file% -C{color_palette} -Z ')  # Create the color palette for the grid file
    os.system(f'gmt grdimage -J{map_projection} -R{region} %grid_file% -Q')  # Plot the data
    os.system(f'gmt psxy {file_name_poly} -W3,red') # Plot the polygon
    os.system(f'gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"{title}" -W0.25p,80/80/80 -Df -V ')  # Plot the coastlines and the title
    os.system(f'gmt text -F+cBL+t"{subtitle}" -N -D6.65c/-1c') # Plot the subtitle
    os.system(f'gmt text -F+cBL+t"{editors}" -N -D5.15c/-1.5c')  # Plot the editors
    # os.system(f'gmt colorbar {colorbar_settings}')  # Plot the colorbar
    if(show_plot):  # Show the plot, if specified
        os.system(f'gmt end show')  # Show the plot
    else:  # Otherwise, just save the plot
        os.system(f'gmt end')  # Save the plot"""

    os.system(f'gmt gmtset FORMAT_GEO_MAP ddd \n set ascii_file={file_name}.txt \n set grid_file={file_name}.grd \n gmt begin {file_name} {img_type} \n gmt xyz2grd %ascii_file% -R{region} -r -I{grid_resolution} -G%grid_file% -V \n gmt grd2cpt %grid_file% -C{color_palette} -Z \n gmt grdimage -J{map_projection} -R{region} %grid_file% -Q \n gmt psxy {file_name_poly} -W3,red \n gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"{title}" -W0.25p,80/80/80 -Df -V \n gmt text -F+cBL+t"{subtitle}" -N -D6.65c/-1c \n gmt text -F+cBL+t"{editors}" -N -D5.15c/-1.5c \n gmt end')
    
    # Move file to the plots folder
    shutil.move(os.path.join(f'{file_name}.{img_type}'), os.path.join('plots', f'{file_name}.{img_type}'))



# Classes
# -----------------------------------------------------------------------------


# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # Plotting of the EWH
    plot_sherical_harmonics(file_name="test_BC_01",
                            img_type="png",
                            grid_resolution="0.5",
                            file_name_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="B-130/65/45/65/18c",
                            region ="-145/-110/45/65",
                            color_palette="haxby",
                            title="British Columbia, Canada",
                            subtitle="Equivalent water heights (EWH)",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V',
                            show_plot=False)

