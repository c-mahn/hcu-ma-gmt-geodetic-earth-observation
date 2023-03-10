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

    command = ""                                                                                                            # Initialize the command
        	
    command += f'gmt gmtset FORMAT_GEO_MAP ddd'                                                                             # Set the format of the map
    command += f' && gmt begin {file_name} {img_type}'                                                                      # Start the plot
    command += f' && gmt xyz2grd ./output/{file_name}.csv -R{region} -r -I{grid_resolution} -G./output/{file_name}.grd -V'  # Convert the ascii file to a grid file
    command += f' && gmt grd2cpt ./output/{file_name}.grd -T-0.25/0.25/0.05 -C{color_palette} -Z'                           # Create the color palette for the grid file
    command += f' && gmt grdimage -J{map_projection} -R{region} ./output/{file_name}.grd -Q'                                # Plot the grid file
    command += f' && gmt psxy {file_poly} -R{region} -W3,red'                                                               # Plot region polygon
    command += f' && gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"{title}" -W0.25p,80/80/80 -Df -N1/1.25p,black -V'                   # Plot the coastline and the title
    command += f' && gmt text -F+cBL+t"{subtitle}" -N -D0c/-1c'                                                             # Plot the subtitle
    command += f' && gmt text -F+cBL+t"{editors}" -N -D0c/-1.5c'                                                            # Plot the editors
    command += f' && gmt colorbar {colorbar_settings}'                                                                      # Plot the colorbar
    command += f' && gmt end'                                                                                               # Save the plot                 

    os.system(command)                                                                                                      # Execute the command
    
    # Move file to the plots folder
    shutil.move(os.path.join(f'{file_name}.{img_type}'), os.path.join('plots', f'{file_name}.{img_type}'))



# Classes
# -----------------------------------------------------------------------------


# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # Plotting of the EWH
    filter_radii = [100, 200, 300, 400, 500]
    """
    for year in range(2003, 2016):
        for month in range(1, 13):
            try:           
                plot_sherical_harmonics(file_name=f"ewh_{year}-{month:02d}_filtered_500km_region_bounding_box",
                                        img_type="png",
                                        grid_resolution="1",
                                        file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                                        map_projection="B-40/70/45/65/18c", # "B-130/65/45/65/18c",
                                        region ="-60/-20/55/85", # "-145/-110/45/65",
                                        color_palette="haxby",
                                        title="groenland", # "British Columbia, Canada",
                                        subtitle=f"Equivalent water heights (EWH) - {year}-{month:02d}, Filterradius:{filter_radius}km",
                                        editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                                        colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                                        )
            except:
                pass
            for filter_radius in filter_radii:
                try:           
                    plot_sherical_harmonics(file_name=f"ewh_{year}-{month:02d}_filtered_{filter_radius}km",
                                            img_type="png",
                                            grid_resolution="1",
                                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                                            map_projection="B-40/70/45/65/18c", # "B-130/65/45/65/18c",
                                            region ="-60/-20/55/85", # "-145/-110/45/65",
                                            color_palette="haxby",
                                            title="groenland", # "British Columbia, Canada",
                                            subtitle=f"Equivalent water heights (EWH) - {year}-{month:02d}, Filterradius:{filter_radius}km",
                                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                                            )
                except:
                    pass
    """



    # British Columbia, Canada


    # Plot the unfiltered EWH
    plot_sherical_harmonics(file_name="ewh_2008-04",
                            img_type="png",
                            grid_resolution="1",
                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="B-130/65/45/65/18c",
                            region ="-145/-110/45/65",
                            color_palette="haxby",
                            title="British Columbia, Canada",
                            subtitle="Unfiltered equivalent water heights (EWH) 2008-04",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                            )

    # Plot the filtered EWH for the different filter radii
    for filter_radius in filter_radii:
        plot_sherical_harmonics(file_name=f"ewh_2008-04_filtered_{filter_radius}km",
                                img_type="png",
                                grid_resolution="1",
                                file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                                map_projection="B-130/65/45/65/18c",
                                region ="-145/-110/45/65",
                                color_palette="haxby",
                                title="British Columbia, Canada",
                                subtitle=f"Gauss-filtered equivalent water heights (EWH) 2008-04, r={filter_radius} km",
                                editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                                colorbar_settings='-Dx0.5c/-2c+w17c/0.35c+h -B0.05+l"EWH [m]" -V'
                                )

    
    """
    # Greenland
    # Plot the unfiltered EWH
    plot_sherical_harmonics(file_name="ewh_2008-04",
                            img_type="png",
                            grid_resolution="1",
                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="B-40/70/45/65/18c",
                            region ="-60/-20/55/85",
                            color_palette="haxby",
                            title="Greenland",
                            subtitle="Unfiltered equivalent water heights (EWH) 2008-04",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                            )
    
    
    # Plot the filtered EWH for the different filter radii
    for filter_radius in filter_radii:
        plot_sherical_harmonics(file_name=f"ewh_2008-04_filtered_{filter_radius}km",
                                img_type="png",
                                grid_resolution="1",
                                file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                                map_projection="B-40/70/45/65/18c",
                                region ="-60/-20/55/85",
                                color_palette="haxby",
                                title="Greenland",
                                subtitle=f"Gauss-filtered equivalent water heights (EWH) 2008-04, r={filter_radius} km",
                                editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                                colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                                )

    
    # World
    # Plot the unfiltered EWH
    plot_sherical_harmonics(file_name="ewh_2008-04",
                            img_type="png",
                            grid_resolution="4",
                            file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                            map_projection="N12c",
                            region ="-180/180/-90/90",
                            color_palette="haxby",
                            title="World",
                            subtitle="Unfiltered equivalent water heights (EWH) 2008-04",
                            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                            )
    
    
    # Plot the filtered EWH for the different filter radii
    for filter_radius in filter_radii:
        plot_sherical_harmonics(file_name=f"ewh_2008-04_filtered_{filter_radius}km",
                                img_type="png",
                                grid_resolution="4",
                                file_poly=os.path.join(main.folder_data, "region_polygon.txt"),
                                map_projection="N12c",
                                region ="-180/180/-90/90",
                                color_palette="haxby",
                                title="World",
                                subtitle=f"Gauss-filtered equivalent water heights (EWH) 2008-04, r={filter_radius} km",
                                editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
                                colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'
                                )
"""