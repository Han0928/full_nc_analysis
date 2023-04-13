import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt

def process_nc_file(filename):
    # Use dask to load the data in chunks
    chunks = {'time': 50, 'model_level_number': 5, 'grid_latitude': 10, 'grid_longitude': 10}
    ds = xr.open_dataset(filename, chunks=chunks)

    # Extract the variable and time data from the dataset
    variable_name = filename.split('/')[-1].split('_')[1:-1]
    variable_data = ds['_'.join(variable_name)]
    variable_data_mean = variable_data.mean(dim=['grid_latitude', 'grid_longitude']).sel(model_level_number=3)
    print(f'variable_data_mean ({variable_name}):', variable_data_mean)
    time_data = ds['time']

    # Compute the value of the mean
    variable_data_mean_value = variable_data_mean.compute()
    time_data_value = time_data.compute()

    # Print the value of the mean
    print(f'The mean of variable_data ({variable_name}) is:', variable_data_mean_value)

    # Convert the time data to a dask array
    time = da.from_array(ds['time'], chunks=chunks['time'])

    return variable_data_mean_value, time_data_value

def process_nc_files(filenames):
    variable_data_mean_values = []
    time_data_values = []

    for filename in filenames:
        variable_data_mean_value, time_data_value = process_nc_file(filename)
        variable_data_mean_values.append(variable_data_mean_value)
        time_data_values.append(time_data_value)

    return variable_data_mean_values, time_data_values

def plot_data(time_data_values, variable_data_mean_values, filenames):
    fig, ax = plt.subplots(figsize=(8,6))
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
    markers = ['o', 's', '^', 'D', 'v']
    
    labels = []
    for filename in filenames:
        variable_name = filename.split('/')[-1].split('_')[8:10]
        variable_name = '_'.join(variable_name)
        label = filename.split('/')[-1].split('_')[3].capitalize() + ' (' + variable_name + ')'
        labels.append(label)

    for i, (time_data_value, variable_data_mean_value) in enumerate(zip(time_data_values, variable_data_mean_values)):
        ax.plot(time_data_value, variable_data_mean_value, label=labels[i], color=colors[i], marker=markers[i], markersize=5, linewidth=1)

    ax.set_xlabel('Time', fontsize=14)
    ax.set_ylabel('# of partiles/air molecule', fontsize=14)
    ax.set_title('Time series of July-August, 2014', fontsize=16, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.legend(labels=labels, fontsize=12, frameon=False)
    plt.xticks(rotation=30)
    plt.show()


# Example usage:
path = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-ct706/full_nc_files/"
filenames = [path+'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
             path+'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc']

variable_data_mean_values, time_data_values = process_nc_files(filenames)

plot_data(time_data_values, variable_data_mean_values, filenames)
