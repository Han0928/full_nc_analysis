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
import datetime
import dateutil.tz
from scipy.interpolate import interp1d, RegularGridInterpolator
import scipy as sp
import netCDF4 as nc

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

# theta to T(k)
def convert_theta_to_temperature(potential_temperature, air_pressure):
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    Rd_cp = 287.05 / 1004.0  
    air_pressure_ratio = air_pressure/p0
    air_pressure_ratio.convert_units('1')
    temperature = potential_temperature*(air_pressure_ratio)**(Rd_cp)  
    return temperature

# from kg/kg to molecule cm-3
def mixing_ratio_to_number_concentration(mixing_ratio_data, air_pressure, actual_temperature):
    zboltz = 1.3807E-23  # (J/K) R = k * N_A, k=J/K, Avogadro's number (N_A)=6.022 x 1023 entities/mol.
    staird = air_pressure / (actual_temperature * zboltz * 1.0E6)  # 1.0E6 from m3 to cm3, another form of ideal gas law
    number_concentration = mixing_ratio_data * staird
    number_concentration.units = 'molecule cm-3'
    return number_concentration

# new version, under testing, for both diameter and number concentration list.
def process_single_file(filename, air_pressure, actual_temperature, bbox, convert_units=True):
    if convert_units:
        variable_name = filename.split('/')[-1].split('_')[1:-1]
        variable_data_cube = iris.load_cube(filename, '_'.join(variable_name))
        variable_data_cube = add_lat_lon(variable_data_cube, bbox)
        number_concentration_data = mixing_ratio_to_number_concentration(variable_data_cube, air_pressure, actual_temperature)
        number_concentration_mean = number_concentration_data.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        number_concentration_mean = number_concentration_mean.extract(iris.Constraint(model_level_number=1))
        
        time_data = variable_data_cube.coord('time')
        # convert time data to list of datetime objects
        time_data_value = time_data.points
        time_units = time_data.units
        epoch = datetime.datetime(1970, 1, 1)
        time_data_datetime = [(epoch + datetime.timedelta(hours=float(tp))) for tp in time_data_value]
        
    else:
        with nc.Dataset(filename, 'r') as ncfile:
            print(ncfile.variables.keys())
            # time_data = ncfile.variables['time']    
        dia_data = iris.load_cube(filename, 'diameter')
        dia_mean = dia_data.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        dia_mean = dia_mean.extract(iris.Constraint(model_level_number=1))
    return number_concentration_mean, time_data_datetime, dia_mean

# new version, under testing, for both diameter and number concentration list.
def process_nc_files(filenames_N, filenames_D, air_pressure, actual_temperature, bbox):
    number_concentration_mean_values = []
    time_data_values = []
    dia_mean_values = []
    # process the nucleation mode files
    for filename in filenames_N:
        number_concentration_mean_value, time_data_value, dia_mean_value = process_single_file(filename, air_pressure, actual_temperature, bbox)
        number_concentration_mean_values.append(number_concentration_mean_value)
        time_data_values.append(time_data_value)
        dia_mean_values.append(dia_mean_value)

    # process the diameter files
    for filename in filenames_D:
        number_concentration_mean_value, time_data_value, dia_mean_value = process_single_file(filename, air_pressure, actual_temperature, bbox, convert_units=False)
        number_concentration_mean_values.append(number_concentration_mean_value)
        # time_data_values.append(time_data_value)
        dia_mean_values.append(dia_mean_value)
    return number_concentration_mean_values, time_data_values, dia_mean_values

#from now on I am testing the Nall's code
def lognormal_cumulative_forcubes(N,r,rbar,sigma):
    total=(N.data/2)*(1+sp.special.erf(np.log(r/rbar.data)/np.sqrt(2)/np.log(sigma)))
    return N.copy(total)


def plot_data(time_data_values, number_concentration_mean_values, filenames_N):
    fig, axes = plt.subplots(5, 1, figsize=(10, 20), sharex=True)
    colors = ['tab:blue', 'tab:orange']
    markers = ['o', 's']
    labels = ['Binary nucleation', 'Updated ion-ternary nucleation']

    for i in range(5):
        axes[i].plot(time_data_values[i], number_concentration_mean_values[i].data, label=labels[0], color=colors[0], marker=markers[0], markersize=3, linewidth=1)
        axes[i].plot(time_data_values[i+5], number_concentration_mean_values[i+5].data, label=labels[1], color=colors[1], marker=markers[1], markersize=3, linewidth=1)
        variable_name = filenames_N[i].split('/')[-1].split('_')[8:10]
        variable_name = '_'.join(variable_name)
        title = filenames_N[i].split('/')[-1].split('_')[3].capitalize() + ' (' + variable_name + ')'
        axes[i].set_title(title, fontsize=14)
        axes[i].tick_params(axis='both', labelsize=12)
        axes[i].set_xlabel('Time', fontsize=12)
        axes[i].set_ylabel('#/cm3', fontsize=12)
        if i == 0:
            axes[i].legend()

    fig.suptitle('BAO tower,40.03(N),105.00(W) \n 0.1 degree, z=1', fontsize=16, fontweight='bold') #https://psl.noaa.gov/technology/bao/site/ 
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.xticks(rotation=30)
    plt.show()
    plt.savefig('output_fig/ncfull_BAO_0.1degree_Z1.png')

path_ct706 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-ct706/full_nc_files/" #i_nuc=2
path_cs093 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-cs093/full_nc_files/" #i_nuc=4

filenames_N = [
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

filenames_D = [
    path_ct706 + 'Rgn_m01s38i401_m01s38i401.nc',
    path_ct706 + 'Rgn_m01s38i402_m01s38i402.nc',
    path_ct706 + 'Rgn_m01s38i403_m01s38i403.nc',
    path_ct706 + 'Rgn_m01s38i404_m01s38i404.nc',
    path_ct706 + 'Rgn_m01s38i405_m01s38i405.nc',
    path_cs093 + 'Rgn_m01s38i405_m01s38i401.nc',
    path_cs093 + 'Rgn_m01s38i405_m01s38i402.nc',
    path_cs093 + 'Rgn_m01s38i405_m01s38i403.nc',
    path_cs093 + 'Rgn_m01s38i405_m01s38i404.nc',
    path_cs093 + 'Rgn_m01s38i405_m01s38i405.nc'
]


potential_temperature_file_ct706 = path_ct706 + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file_ct706 = path_ct706 + 'Rgn_air_pressure_m01s00i408.nc'

potential_temperature_file_cs093 = path_cs093 + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file_cs093 = path_cs093 + 'Rgn_air_pressure_m01s00i408.nc'

bbox = [-105.05, -104.95, 39.95, 40.05] #BAO tower,40 03 00.10028(N) Longitude: 105 00 13.80781(W) Elevation: 1584 m Height: 985 ft (300 m)
# bbox = [-106.65, -106.55, 40.40, 40.50] #Storm Peak Lab:40.45° N, 106.6° W, high altitude
potential_temperature_ct706, air_pressure_ct706 = read_pt_data(potential_temperature_file_ct706, air_pressure_file_ct706, bbox)
actual_temperature_ct706 = convert_theta_to_temperature(potential_temperature_ct706, air_pressure_ct706)

potential_temperature_cs093, air_pressure_cs093 = read_pt_data(potential_temperature_file_cs093, air_pressure_file_cs093, bbox)
actual_temperature_cs093 = convert_theta_to_temperature(potential_temperature_cs093, air_pressure_cs093)

number_concentration_mean_values_ct706, time_data_values_ct706,dia_mean_values_ct706 = process_nc_files(filenames_N[:5], filenames_D[:5],air_pressure_ct706, actual_temperature_ct706, bbox)
number_concentration_mean_values_cs093, time_data_values_cs093,dia_mean_values_cs093 = process_nc_files(filenames_N[5:], filenames_D[5:],air_pressure_cs093, actual_temperature_cs093, bbox)

number_concentration_mean_values = number_concentration_mean_values_ct706 + number_concentration_mean_values_cs093
time_data_values = time_data_values_ct706 + time_data_values_cs093

plot_data(time_data_values, number_concentration_mean_values, filenames_N[:5] + filenames_N[5:])

