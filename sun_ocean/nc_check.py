import netCDF4 as nc
import numpy as np

# Open the netCDF file
nc_file = nc.Dataset('pacific_ocean_large.nc')

# Access the variable of interest
u_var = nc_file.variables['u']
print('u_var', u_var)
# Get variable attributes
units = u_var.units
fill_value = u_var._FillValue  # Assuming the fill value is stored in the _FillValue attribute
print(f"Units: {units}")
print(f"Fill value: {fill_value}")
# Get the raw data
u_data = u_var[:]

# Find indices of large values in the range of -32xxx
large_value_indices = np.where((u_data >= -32000) & (u_data < -31999))

# Print information about the large values
for index in np.transpose(large_value_indices):
    time_index, pressure_index, lat_index, lon_index = index
    value = u_data[time_index, pressure_index, lat_index, lon_index]
    print(f"Large value at time index {time_index}, pressure index {pressure_index}, latitude index {lat_index}, longitude index {lon_index}: {value}")

# Close the netCDF file
nc_file.close()
