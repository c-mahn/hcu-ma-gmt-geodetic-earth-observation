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
                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
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

    command = ""
    command += f'gmt gmtset FORMAT_GEO_MAP ddd'     # Set the format of the map
    command += f' && gmt begin {file_name} {img_type}'  # Start the plot
    os.system(f'gmt gmtset FORMAT_GEO_MAP ddd')     # Set the format of the map
    os.system(f'gmt begin {file_name} {img_type}')  # Start the plot
    os.system(f'gmt xyz2grd ./output/{file_name}.csv -R{region} -r -I{grid_resolution} -G -V')           # Convert the ascii file to a grid file
    os.system(f'gmt grd2cpt ./output/{file_name}.grd -C{color_palette} -Z')                                          # Create the color palette for the grid file
    os.system(f'gmt grdimage -J{map_projection} -R{region} ./output/{file_name}.grd -Q')                             # Plot the grid file                    
    os.system(f'gmt psxy {file_poly} -R{region} -W3,red')                                                          # Plot region polygon
    os.system(f'gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"{title}" -W0.25p,80/80/80 -Df -N1/1.25p,black -V')   # Plot the coastline and the title
    os.system(f'gmt text -F+cBL+t"{subtitle}" -N -D6.65c/-1c')  # Plot the subtitle
    os.system(f'gmt text -F+cBL+t"{editors}" -N -D5.15c/-1.5c') # Plot the editors
    os.system(f'gmt colorbar {colorbar_settings}')               # Plot the colorbar
    os.system(f'gmt end')
        
    if(show_plot):  # Show the plot, if specified
        os.system(f'gmt end show')  # Show the plot
    else:  # Otherwise, just save the plot
        os.system(f'gmt end')  # Save the plot
    
    # Move file to the plots folder
    shutil.move(os.path.join(f'{file_name}.{img_type}'), os.path.join('plots', f'{file_name}.{img_type}'))



# Classes
# -----------------------------------------------------------------------------


# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # Plotting of the EWH
    plot_sherical_harmonics(file_name="monthly_ewh_filtered_300km_2014-09",
                            img_type="png",
                            grid_resolution="1",
                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="B-130/65/45/65/18c",
                            region ="-145/-110/45/65",
                            color_palette="haxby",
                            title="British Columbia, Canada",
                            subtitle="Equivalent water heights (EWH)",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V',
                            show_plot=False)
"""
    plot_sherical_harmonics(file_name="test_BC_01",
                            img_type="png",
                            grid_resolution="0.5",
                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="B-130/65/45/65/18c",
                            region ="-145/-110/45/65",
                            color_palette="haxby",
                            title="British Columbia, Canada",
                            subtitle="Equivalent water heights (EWH)",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V',
                            show_plot=False)
"""