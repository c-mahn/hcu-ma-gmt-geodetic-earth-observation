# Geodetic Earth Observation - Task 03
# #############################################################################

# Authors:
# Christopher Mahn
# Silas Teske
# Joshua Wolf

# #############################################################################

# Import of Libraries
# -----------------------------------------------------------------------------

import numpy as np
import os
import shutil
import netCDF4 as nc
import functions as fu

# #############################################################################


# Constants
# -----------------------------------------------------------------------------

g                   = 9.807                 # Gravitational acceleration [m/s^2]
omega               = (2*np.pi)/(24*60*60)  # Rotational velocity of the Earth
radius              = 6.378137000e+06       # m
rho_grad            = 180/np.pi             # Conversion factor from radians to degrees

# Functions
# -----------------------------------------------------------------------------

def run_gmt(input_file_name="test_input",
            output_file_name=None,
            sample_and_cut=False,
            print_vector=False,
            vector_1=None,
            vector_2=None,
            img_type="png",
            grid_resolution="2m",
            map_projection="B-80/40/20/40/18c",         # Florida ("B-130/65/45/65/18c" British Columbia)
            region = "-90/-70/20/40",                   # Florida ("-145/-110/45/65" British Columbia)
            color_palette="haxby",
            color_settings="-1/1/0.001",
            title="Florida, USA",                       # Florida ("British Columbia, Canada" British Columbia)
            subtitle=None,
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V'):
    """
    This function runs GMT commands to sample and cut a grid, plot a grid, and plot velocity vectors.

    Args:
    input_file_name (str): Name of the input file
    output_file_name (str): Name of the output file
    sample_and_cut (bool): If True, the grid will be sampled and cut
    print_vector (bool): If True, the velocity vectors will be plotted
    vector_1 (str): Name of the first vector file
    vector_2 (str): Name of the second vector file
    img_type (str): Type of image to be plotted
    grid_resolution (str): Resolution of the grid
    map_projection (str): Map projection
    region (str): Region of the map or the grid to be cut
    color_palette (str): Color palette
    color_settings (str): Color settings
    title (str): Title of the plot
    subtitle (str): Subtitle of the plot
    editors (str): Name of the editors
    colorbar_settings (str): Colorbar settings

    """
    if output_file_name is None:
        output_file_name = input_file_name
    
    if subtitle is None:
        subtitle = output_file_name

    command = ""
    command += f'gmt gmtset FORMAT_GEO_MAP ddd'                                                 # Set the format of the map
    if(print_vector is True):
        command += f' && gmt begin vectors_{output_file_name} {img_type}'                       # Start the plot
    elif(sample_and_cut is True):
        command += f' && gmt begin'                                                             # Dont start the plot
    else:
        command += f' && gmt begin {output_file_name} {img_type}'                               # Start the plot
    if(sample_and_cut is True):
        command += f' && gmt grdsample ./data/{input_file_name}.nc -G./data/{input_file_name}_{grid_resolution}.nc -I{grid_resolution}'     # Sample the grid
        command += f' && gmt grdcut ./data/{input_file_name}_{grid_resolution}.nc -G./data/{output_file_name}.nc -R{region}'                # Cut the grid
    command += f' && gmt grdinfo ./data/{output_file_name}.nc'                                                                              # Get the information of the grid
    command += f' && gmt grd2cpt ./data/{output_file_name}.nc -T{color_settings} -C{color_palette} -Z'                                      # Apply the color palette
    if(print_vector is True):
        command += f' && gmt grdimage -J{map_projection} -R{region} ./data/{output_file_name}.nc -Q'                    # Plot the grid
        command += f' && gmt grdvector ./data/{vector_1}.nc ./data/{vector_2}.nc -WRed -Ix2 -S2 -Q0.2+e'                # Plot the vectors
    else:
        command += f' && gmt grdimage -J{map_projection} -R{region} ./data/{output_file_name}.nc -Q'                    # Plot the grid
    command += f' && gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"{title}" -W0.25p,80/80/80 -Df -N1/1.25p,black -V '              # Plot the coastlines and the title
    command += f' && gmt text -F+cBL+t"{subtitle}" -N -D0c/-1c'                                                         # Plot the subtitle
    command += f' && gmt text -F+cBL+t"{editors}" -N -D0c/-1.5c'                                                        # Plot the editors
    command += f' && gmt colorbar {colorbar_settings}'                                                                  # Plot the colorbar
    command += f' && gmt end'                                                                                           # Save the plot                 

    os.system(command)
    
    # Move image file to the plots folder
    if(print_vector is True):
        shutil.move(os.path.join(f'vectors_{output_file_name}.{img_type}'), os.path.join('plots', f'vectors_{output_file_name}.{img_type}'))
    elif(sample_and_cut is True):
        pass
    else:
        shutil.move(os.path.join(f'{output_file_name}.{img_type}'), os.path.join('plots', f'{output_file_name}.{img_type}'))


def calc_velocity(mdt, lon, lat, grid_spacing=0.5):
    u = np.zeros((len(lat), len(lon)))
    v = np.zeros((len(lat), len(lon)))
    grid_spacing /= rho_grad
    for index_phi, phi in enumerate(lat):
        phi /= rho_grad
        f = 2*omega*np.sin(phi)
        for index_lam, lam in enumerate(lon):
            lam /= rho_grad
            if index_lam == 0:
                continue
            elif index_lam == len(lon)-1:
                continue
            elif index_phi == 0:
                continue
            elif index_phi == len(lat)-1:
                continue
            else:
                dx = radius * np.cos(phi)*grid_spacing
                dy = radius * grid_spacing
                dhx = mdt[index_phi][index_lam + 1] - mdt[index_phi][index_lam - 1]
                dhy = mdt[index_phi + 1][index_lam] - mdt[index_phi - 1][index_lam]
                u[index_phi][index_lam] = -(g/f)*(dhy/(2*dy))
                v[index_phi][index_lam] = (g/f)*(dhx/(2*dx))
    return u, v

# Beginning of the programm
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Resampling and cutting the grid-data

    run_gmt(input_file_name="DTU15MDT_2min.mdt",
            output_file_name="DTU15MDT_30min_cut.mdt",
            sample_and_cut=True,
            grid_resolution="30m",
            )
 
    # Plotting the grid-data saved with the function run_gmt

    run_gmt(input_file_name="DTU15MDT_30min_cut.mdt",
            grid_resolution="30m",
            color_settings="-1/1/0.2",
            subtitle="Magnitude of geostrophic surface velocity",
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.5+l"MDT [m]" -V'
            )
     
    # Importing the grid-data with netCDF4 and saving it with the function save_grid

    data_mdt = nc.Dataset(os.path.join(f'./data/DTU15MDT_30min_cut.mdt.nc'))
    lon = data_mdt.variables["lon"][:]
    lat = data_mdt.variables["lat"][:]
    mdt = data_mdt.variables["z"][:]
    fu.save_grid(os.path.join(f'./data/DTU15MDT_30min_cut_pyout.mdt.nc'), mdt, lon, lat)

    # Calculation of geostrophic currents

    u, v = calc_velocity(mdt, lon, lat, grid_spacing=0.5)

    # Saving the geostrophic currents

    fu.save_grid(os.path.join(f'./data/DTU15MDT_30min_velocity.mdt.nc'), np.sqrt(u**2+v**2), lon, lat)

    fu.save_grid(os.path.join(f'./data/velocity_components_east.mdt.nc'), u, lon, lat)
    fu.save_grid(os.path.join(f'./data/velocity_components_north.mdt.nc'), v, lon, lat)

    # Plotting the geostrophic current vectors

    run_gmt(input_file_name="DTU15MDT_30min_velocity.mdt",
            print_vector=True,
            vector_1="velocity_components_east.mdt",
            vector_2="velocity_components_north.mdt",
            grid_resolution="30m",
            color_settings="0/0.5/0.1",
            subtitle="Ocean current velocities with geostrophic velocity vectors",
            editors="Editors: Christopher Mahn, Silas Teske, Joshua Wolf",
            colorbar_settings='-Dx0c/-2c+w17c/0.35c+h -B0.2+l"Velocity [m/s]" -V'
            )
    