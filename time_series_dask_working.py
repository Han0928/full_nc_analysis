import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import os
import matplotlib.dates as mdates

def read_pt_data(potential_temperature_file, air_pressure_file):
    chunks = {'time': 20, 'model_level_number': 5, 'grid_latitude': 10, 'grid_longitude': 10}
    
    ds_theta = xr.open_dataset(potential_temperature_file, chunks=chunks)
    ds_p = xr.open_dataset(air_pressure_file, chunks=chunks)

    potential_temperature = ds_theta['air_potential_temperature']
    air_pressure = ds_p['air_pressure']
    
    return potential_temperature, air_pressure

def convert_theta_to_temperature(potential_temperature, air_pressure):
    Rd_cp = 287.05 / 1004.0
    temperature = potential_temperature * (air_pressure / 100000.0) ** Rd_cp  #hpa? 
    T = temperature.mean().compute()
    print('T:', T) 

    return temperature


def mixing_ratio_to_number_concentration(mixing_ratio_data, air_pressure, actual_temperature):
#     R = 287.0
#     Na = 6.022e23
#     air_density = air_pressure / (R * actual_temperature)
#     number_concentration = mixing_ratio_data * air_density * Na
#     return number_concentration

# Now I want to test the unit conversion, the above commented code works(and give me #/m3)
    zboltz=1.3807E-23
    staird=air_pressure/(actual_temperature*zboltz*1.0E6) # 1.0E6 from m3 to cm3
    number_concentration = mixing_ratio_data * staird
    return number_concentration   

def process_nc_file(filename, air_pressure, actual_temperature):
    chunks = {'time': 20, 'model_level_number': 5, 'grid_latitude': 10, 'grid_longitude': 10}
    ds = xr.open_dataset(filename, chunks=chunks)
    variable_name = filename.split('/')[-1].split('_')[1:-1]
    variable_data = ds['_'.join(variable_name)]
    number_concentration_data = mixing_ratio_to_number_concentration(variable_data, air_pressure, actual_temperature)
    number_concentration_mean = number_concentration_data.mean(dim=['grid_latitude', 'grid_longitude']).sel(model_level_number=2)
    time_data = ds['time']
    number_concentration_mean_value = number_concentration_mean.compute()
    time_data_value = time_data.compute()
    return number_concentration_mean_value, time_data_value

#3 is too high, maybe 2 is better, given the tower height.

def process_nc_files(filenames, air_pressure, actual_temperature):
    number_concentration_mean_values = []
    time_data_values = []
    for filename in filenames:
        number_concentration_mean_value, time_data_value = process_nc_file(filename, air_pressure, actual_temperature)
        number_concentration_mean_values.append(number_concentration_mean_value)
        time_data_values.append(time_data_value)
    return number_concentration_mean_values, time_data_values

# Updated plot_data function
def plot_data(time_data_values, number_concentration_mean_values, filenames):
    fig, axes = plt.subplots(5, 1, figsize=(6, 20), sharex=True)
    colors = ['tab:blue', 'tab:orange']
    markers = ['o', 's']

    labels = ['Binary nucleation', 'Updated ion-ternary nucleation']

    for i in range(5):
        axes[i].plot(time_data_values[i], number_concentration_mean_values[i], label=labels[0], color=colors[0], marker=markers[0], markersize=3, linewidth=1)
        axes[i].plot(time_data_values[i+5], number_concentration_mean_values[i+5], label=labels[1], color=colors[1], marker=markers[1], markersize=3, linewidth=1)

        variable_name = filenames[i].split('/')[-1].split('_')[8:10]
        variable_name = '_'.join(variable_name)
        title = filenames[i].split('/')[-1].split('_')[3].capitalize() + ' (' + variable_name + ')'
        axes[i].set_title(title, fontsize=14)

        axes[i].tick_params(axis='both', labelsize=12)
        axes[i].set_xlabel('Time', fontsize=12)
        axes[i].set_ylabel('#/cm3', fontsize=12)
        axes[i].legend()

        # Set date labels every 2 days
        axes[i].xaxis.set_major_locator(mdates.DayLocator(interval=2))
        axes[i].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    fig.suptitle('Number Conc, July-August, 2014', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.xticks(rotation=30)
    plt.show()

# Now need to read in the file
path_ct706 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-ct706/full_nc_files/" #i_nuc=2
path_cs093 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-cs093/full_nc_files/" #i_nuc=4

filenames = [
    # u-ct706 files
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc',
    # u-cs093 files
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

potential_temperature_ct706, air_pressure_ct706 = read_pt_data(potential_temperature_file_ct706, air_pressure_file_ct706)
actual_temperature_ct706 = convert_theta_to_temperature(potential_temperature_ct706, air_pressure_ct706)

potential_temperature_cs093, air_pressure_cs093 = read_pt_data(potential_temperature_file_cs093, air_pressure_file_cs093)
actual_temperature_cs093 = convert_theta_to_temperature(potential_temperature_cs093, air_pressure_cs093)

number_concentration_mean_values_ct706, time_data_values_ct706 = process_nc_files(filenames[:5], air_pressure_ct706, actual_temperature_ct706)
number_concentration_mean_values_cs093, time_data_values_cs093 = process_nc_files(filenames[5:], air_pressure_cs093, actual_temperature_cs093)

number_concentration_mean_values = number_concentration_mean_values_ct706 + number_concentration_mean_values_cs093
time_data_values = time_data_values_ct706 + time_data_values_cs093

plot_data(time_data_values, number_concentration_mean_values, filenames[:5] + filenames[5:])


