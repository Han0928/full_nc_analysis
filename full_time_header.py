import xarray as xr

filename = 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc'
ds = xr.open_dataset(filename)
print(ds['time'].attrs)

# Extract the "time" variable and convert to a pandas datetime object
time = ds['time'].values
#time_units = ds['time'].units
time_units = ds['time'].attrs['units']

time_calendar = ds['time'].calendar
pandas_time = xr.coding.times.decode_cf_datetime(time, units=time_units, calendar=time_calendar)

# Convert the pandas datetime object to a human-readable format
human_readable_time = pandas_time.strftime('%Y-%m-%d %H:%M:%S')

# Print the human-readable time
print(f'The time coordinate of the NetCDF file is: {human_readable_time}')

