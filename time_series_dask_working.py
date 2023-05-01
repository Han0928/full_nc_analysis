# %load /ocean/projects/atm200005p/ding0928/script_full_nc/time_series_dask_working.py
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import os
import matplotlib.dates as mdates
import numpy as np
import iris
import iris.coord_systems as cs
import iris.coord_systems as coord_systems

os.environ["OPENBLAS_NUM_THREADS"] = "8"

def bbox_extract_2Dcoords(cube, bbox):
    minmax = lambda x: (np.min(x), np.max(x))
    lons = cube.coord('longitude').points 
    lats = cube.coord('latitude').points
    inregion = np.logical_and(np.logical_and(lons > bbox[0],
                                             lons < bbox[1]),
                              np.logical_and(lats > bbox[2],
                                             lats < bbox[3]))
    region_inds = np.where(inregion)
    imin, imax = minmax(region_inds[0])
    jmin, jmax = minmax(region_inds[1])
    return cube[..., imin:imax+1, jmin:jmax+1]

def add_lat_lon(cube, bbox):
    polelat = cube.coord('grid_longitude').coord_system.grid_north_pole_latitude
    polelon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude  

    source_lon = cube.coord('grid_longitude').points
    source_lat = cube.coord('grid_latitude').points
    lat2d = np.transpose(np.tile(source_lat,[len(source_lon),1]))
    lon2d = np.tile(source_lon,[len(source_lat),1])

    lons, lats = iris.analysis.cartography.unrotate_pole(lon2d, lat2d, polelon, polelat)

    longit = iris.coords.AuxCoord(lons,'longitude', units='degrees', coord_system=cs.GeogCS(6371229.0))
    latit =  iris.coords.AuxCoord(lats,'latitude', units='degrees', coord_system=cs.GeogCS(6371229.0))
    
    i_test = 1 # a label for turn-on/off 
    if i_test == 0:
        # Determine the dimensions to add the latitude and longitude coordinates
        dims = tuple(range(cube.ndim))
        for dim in ('time', 'model_level_number', 'grid_latitude', 'grid_longitude'):
            if dim in cube.dim_coords:
                dims = tuple(d for d in dims if d != cube.coord_dims(dim)[0])

    cube.add_aux_coord(longit, (2,3)) 
    cube.add_aux_coord(latit, (2,3))

    return bbox_extract_2Dcoords(cube, bbox)

def read_pt_data(potential_temperature_file, air_pressure_file, bbox):
    potential_temperature_cube = iris.load_cube(potential_temperature_file)
    air_pressure_cube = iris.load_cube(air_pressure_file)
    print(potential_temperature_cube.coord('grid_longitude').points.min(), potential_temperature_cube.coord('grid_longitude').points.max())
    print(potential_temperature_cube.coord('grid_latitude').points.min(), potential_temperature_cube.coord('grid_latitude').points.max())
    print(potential_temperature_cube.units) # K
    print(air_pressure_cube.units) #Pa

    # Add the latitude and longitude coordinates to the cubes
    potential_temperature_cube = add_lat_lon(potential_temperature_cube, bbox)
    air_pressure_cube = add_lat_lon(air_pressure_cube, bbox)
    return potential_temperature_cube, air_pressure_cube

# a subroutine to convert theta to T(k)
def convert_theta_to_temperature(potential_temperature, air_pressure):
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    Rd_cp = 287.05 / 1004.0  
    air_pressure_ratio = air_pressure/p0
    air_pressure_ratio.convert_units('1')
    temperature = potential_temperature*(air_pressure_ratio)**(Rd_cp)
    # now the T looks correct, in 280~, so I commented it out now
#     for i, temp in enumerate(temperature.data.flatten()):
#         print(f'Temperature at grid point {i+1}: {temp:.2f} K')       
    return temperature

# to convert from kg/kg to molecule cm-3
def mixing_ratio_to_number_concentration(mixing_ratio_data, air_pressure, actual_temperature):
    zboltz = 1.3807E-23  # (J/K) R = k * N_A, k=J/K, Avogadro's number (N_A)=6.022 x 1023 entities/mol.
    staird = air_pressure / (actual_temperature * zboltz * 1.0E6)  # 1.0E6 from m3 to cm3, another form of ideal gas law
    number_concentration = mixing_ratio_data * staird
    number_concentration.units = 'molecule cm-3'
    return number_concentration

#iris to process one single file
def process_single_file(filename, air_pressure, actual_temperature, bbox):
    variable_name = filename.split('/')[-1].split('_')[1:-1]
    variable_data_cube = iris.load_cube(filename, '_'.join(variable_name))
    variable_data_cube = add_lat_lon(variable_data_cube, bbox)

    number_concentration_data = mixing_ratio_to_number_concentration(variable_data_cube, air_pressure, actual_temperature)
    number_concentration_mean = number_concentration_data.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
    number_concentration_mean = number_concentration_mean.extract(iris.Constraint(model_level_number=2))

    time_data = variable_data_cube.coord('time')
    time_data_value = time_data.points
    return number_concentration_mean, time_data_value

# to process all the files within a loop
def process_nc_files(filenames, air_pressure, actual_temperature, bbox):
    number_concentration_mean_values = []
    time_data_values = []
    for filename in filenames:
        number_concentration_mean_value, time_data_value = process_single_file(filename, air_pressure, actual_temperature, bbox)
        number_concentration_mean_values.append(number_concentration_mean_value)
        time_data_values.append(time_data_value)
    return number_concentration_mean_values, time_data_values


def plot_data(time_data_values, number_concentration_mean_values, filenames):
    fig, axes = plt.subplots(5, 1, figsize=(6, 20), sharex=True)
    colors = ['tab:blue', 'tab:orange']
    markers = ['o', 's']
    labels = ['Binary nucleation', 'Updated ion-ternary nucleation']
    
    for i in range(5):
        # print("time_data_values[i].dimension",time_data_values[i])
        # print("number_concentration_mean_values[i][0].data.dimension",number_concentration_mean_values[i].data)
        axes[i].plot(time_data_values[i], number_concentration_mean_values[i].data, label=labels[0], color=colors[0], marker=markers[0], markersize=3, linewidth=1)
        axes[i].plot(time_data_values[i+5], number_concentration_mean_values[i+5].data, label=labels[1], color=colors[1], marker=markers[1], markersize=3, linewidth=1)

        variable_name = filenames[i].split('/')[-1].split('_')[8:10]
        variable_name = '_'.join(variable_name)
        title = filenames[i].split('/')[-1].split('_')[3].capitalize() + ' (' + variable_name + ')'
        axes[i].set_title(title, fontsize=14)

        axes[i].tick_params(axis='both', labelsize=12)
        axes[i].set_xlabel('Time', fontsize=12)
        axes[i].set_ylabel('#/cm3', fontsize=12)
        axes[i].legend()

    fig.suptitle('BAO tower:40.5°N, -105W', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.xticks(rotation=30)
    plt.show()
    plt.savefig('BAO_tower_40.5°N_-105.png')

# Now need to read in the file
path_ct706 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-ct706/full_nc_files/" #i_nuc=2
path_cs093 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-cs093/full_nc_files/" #i_nuc=4

filenames = [
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc'
]

potential_temperature_file_ct706 = path_ct706 + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file_ct706 = path_ct706 + 'Rgn_air_pressure_m01s00i408.nc'

potential_temperature_file_cs093 = path_cs093 + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file_cs093 = path_cs093 + 'Rgn_air_pressure_m01s00i408.nc'

# Define the bounding box (in degrees) for the area of interest
bbox = [-105, -104.5, 40.4, 40.8] #BAO tower, low altitude
#bbox = [-107, -106.5, 40.5, 40.8] #Storm Peak Lab:40.45° N, 106.6°  wester-high-altitude
potential_temperature_ct706, air_pressure_ct706 = read_pt_data(potential_temperature_file_ct706, air_pressure_file_ct706, bbox)
actual_temperature_ct706 = convert_theta_to_temperature(potential_temperature_ct706, air_pressure_ct706)

potential_temperature_cs093, air_pressure_cs093 = read_pt_data(potential_temperature_file_cs093, air_pressure_file_cs093, bbox)
actual_temperature_cs093 = convert_theta_to_temperature(potential_temperature_cs093, air_pressure_cs093)

number_concentration_mean_values_ct706, time_data_values_ct706 = process_nc_files(filenames[:5], air_pressure_ct706, actual_temperature_ct706, bbox)
number_concentration_mean_values_cs093, time_data_values_cs093 = process_nc_files(filenames[5:], air_pressure_cs093, actual_temperature_cs093, bbox)
number_concentration_mean_values = number_concentration_mean_values_ct706 + number_concentration_mean_values_cs093
time_data_values = time_data_values_ct706 + time_data_values_cs093

plot_data(time_data_values, number_concentration_mean_values, filenames[:5] + filenames[5:])


