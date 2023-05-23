import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np

# Open the netCDF file
nc_file = nc.Dataset('pacific_ocean.nc')

# Read the necessary variables
u = nc_file.variables['u'][:]
v = nc_file.variables['v'][:]
latitude = nc_file.variables['latitude'][:]
longitude = nc_file.variables['longitude'][:]

# Compute wind speed from u and v components
wind_speed = np.sqrt(u**2 + v**2)

# Set the domain for the Pacific Ocean area
lon_min, lon_max = -180, -100
lat_min, lat_max = 0, 60

# Find the indices corresponding to the desired domain
lon_indices = np.where((longitude >= lon_min) & (longitude <= lon_max))[0]
lat_indices = np.where((latitude >= lat_min) & (latitude <= lat_max))[0]

# Select the wind speed data within the desired domain
wind_speed_domain = wind_speed[:, :, :, lat_indices, lon_indices]

# Plotting the spatial map of wind speed
lon_grid, lat_grid = np.meshgrid(longitude[lon_indices], latitude[lat_indices])

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot wind speed as a filled contour map
c = ax.contourf(lon_grid, lat_grid, wind_speed_domain[0, 0, 0], cmap='coolwarm')
plt.colorbar(c, label='Wind Speed (m/s)')

# Add a marker for the specific point (9°N, 156°W)
ax.plot(-156, 9, 'ro', markersize=5, transform=ccrs.PlateCarree())

# Add coastlines and gridlines
ax.coastlines()
ax.gridlines(draw_labels=True)

# Set plot title and labels
plt.title('Spatial Map of Wind Speed in the Pacific Ocean')
plt.xlabel('Longitude (degrees_east)')
plt.ylabel('Latitude (degrees_north)')

plt.show()

# Close the netCDF file
nc_file.close()
