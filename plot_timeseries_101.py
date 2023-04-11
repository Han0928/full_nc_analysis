import xarray as xr
import matplotlib.pyplot as plt

# Open the NetCDF file and load the data
filename = 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc'
ds = xr.open_dataset(filename)

# Extract the variable and time data from the dataset
variable_data = ds['number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air'].values
time_data = ds['time'].values

# Convert the time data to a pandas datetime object
time_units = ds['time'].attrs['units']
time_calendar = ds['time'].calendar
time = xr.coding.times.decode_cf_datetime(time_data, units=time_units, calendar=time_calendar)

# Create a linear plot of the variable against time
plt.plot(time, variable_data)
plt.xlabel('Time')
plt.ylabel('Variable')
plt.title('Linear Plot of Variable against Time')
plt.show()

