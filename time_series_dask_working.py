import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import os

def read_pt_data(potential_temperature_file, air_pressure_file):
    chunks = {'time': 20, 'model_level_number': 5, 'grid_latitude': 10, 'grid_longitude': 10}
    
    ds_theta = xr.open_dataset(potential_temperature_file, chunks=chunks)
    ds_p = xr.open_dataset(air_pressure_file, chunks=chunks)

    potential_temperature = ds_theta['air_potential_temperature']
    air_pressure = ds_p['air_pressure']
    
    return potential_temperature, air_pressure

def convert_theta_to_temperature(potential_temperature, air_pressure):
    Rd_cp = 287.05 / 1004.0
    temperature = potential_temperature * (air_pressure / 100000.0) ** Rd_cp
    return temperature

def mixing_ratio_to_number_concentration(mixing_ratio_data, air_pressure, actual_temperature):
    R = 287.0
    Na = 6.022e23
    air_density = air_pressure / (R * actual_temperature)
    number_concentration = mixing_ratio_data * air_density * Na
    return number_concentration

def process_nc_file(filename, air_pressure, actual_temperature):
    chunks = {'time': 20, 'model_level_number': 5, 'grid_latitude': 10, 'grid_longitude': 10}
    ds = xr.open_dataset(filename, chunks=chunks)
    variable_name = filename.split('/')[-1].split('_')[1:-1]
    variable_data = ds['_'.join(variable_name)]
    number_concentration_data = mixing_ratio_to_number_concentration(variable_data, air_pressure, actual_temperature)
    number_concentration_mean = number_concentration_data.mean(dim=['grid_latitude', 'grid_longitude']).sel(model_level_number=3)
    time_data = ds['time']
    number_concentration_mean_value = number_concentration_mean.compute()
    time_data_value = time_data.compute()
    return number_concentration_mean_value, time_data_value

def process_nc_files(filenames, air_pressure, actual_temperature):
    number_concentration_mean_values = []
    time_data_values = []
    for filename in filenames:
        number_concentration_mean_value, time_data_value = process_nc_file(filename, air_pressure, actual_temperature)
        number_concentration_mean_values.append(number_concentration_mean_value)
        time_data_values.append(time_data_value)
    return number_concentration_mean_values, time_data_values

def plot_data(time_data_values, number_concentration_mean_values, filenames):
    fig, ax = plt.subplots()
    for i, (time_data, number_concentration_mean) in enumerate(zip(time_data_values, number_concentration_mean_values)):
        ax.plot(time_data, number_concentration_mean, label=filenames[i].split('/')[-1])
    ax.legend()
    plt.xlabel('Time')
    plt.ylabel('Number concentration mean')
    plt.show()

# Example usage:
path = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-ct706/full_nc_files/"
filenames = [path+'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc']

potential_temperature_file = path + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file = path + 'Rgn_air_pressure_m01s00i408.nc'

potential_temperature, air_pressure = read_pt_data(potential_temperature_file, air_pressure_file)
actual_temperature = convert_theta_to_temperature(potential_temperature, air_pressure)

number_concentration_mean_values, time_data_values = process_nc_files(filenames, air_pressure, actual_temperature)
plot_data(time_data_values, number_concentration_mean_values, filenames)
