import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from PIL import Image
import calendar

def plot_spatial_map_direction(u, v, latitude, longitude, month):
    wind_speed = np.sqrt(u**2 + v**2)
    wind_direction = np.arctan2(v, u) * 180 / np.pi

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-160, -150, 6, 12], ccrs.PlateCarree())
    ax.coastlines(linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.3)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    # Draw the region of interest box
    lon_min, lat_min = -157, 8
    lon_max, lat_max = -155, 10
    ax.plot([lon_min, lon_min, lon_max, lon_max, lon_min], [lat_min, lat_max, lat_max, lat_min, lat_min],
            transform=ccrs.PlateCarree(), color='red', linewidth=2)

    label = '*Station\n(-156W, 9N))'
    ax.text(-156, 9, label, fontsize=12, color='black', ha='center', va='center', transform=ccrs.PlateCarree())

    ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    ax.gridlines(draw_labels=False, linewidth=0.0, linestyle='--')

    # Plot wind speed direction with thinner arrow thickness
    im = ax.quiver(longitude, latitude, u, v, wind_speed, cmap='jet', transform=ccrs.PlateCarree(),
                   linewidth=0.2, clim=[np.min(wind_speed), np.max(wind_speed)])
    # Set the arrow angles to match wind direction
    im.set_UVC(wind_speed * np.cos(wind_direction), wind_speed * np.sin(wind_direction))

    # Remove right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', shrink=0.7)
    cbar.set_label('Wind Speed (m/s)', fontsize=12)
    plt.title(f'Wind Direction\nMonth: {calendar.month_name[month]}', fontsize=12)

# Open the netCDF file
nc_file = nc.Dataset('pacific_ocean_large.nc')
latitude = nc_file.variables['latitude'][:]
longitude = nc_file.variables['longitude'][:]
u = nc_file.variables['u'][:]
v = nc_file.variables['v'][:]
times = nc_file.variables['time'][:]  # Get all time slots

# Convert time to a human-readable format
base_date = datetime.datetime(1900, 1, 1)
converted_times = [base_date + datetime.timedelta(hours=int(time)) for time in times]
# Create a directory to store the images
os.makedirs('spatial_maps', exist_ok=True)

# Calculate the monthly averages
monthly_accumulated_u = np.zeros((12, u.shape[1], u.shape[2], u.shape[3], u.shape[4]))
monthly_accumulated_v = np.zeros((12, u.shape[1], u.shape[2], u.shape[3], u.shape[4]))
monthly_counts = np.zeros((12,), dtype=int)

for i, time in enumerate(converted_times):
    # Filter out data from 2023
    if time.year == 2023:
        continue
    
    month = time.month - 1  # Convert to 0-based index
    print('converted_times:', converted_times)
    print('time.month:', time.month)
    print("Before accumulation - Month:", month+1)
    print("----------------------")
    month = time.month - 1  # Convert to 0-based index
    monthly_accumulated_u[month] += u[i]
    monthly_accumulated_v[month] += v[i]
    monthly_counts[month] += 1
    print("After accumulation - Month:", month+1)
    print("Counts:")
    print(monthly_counts[month])
    print("----------------------")

# Calculate the monthly averages by dividing by the counts
monthly_average_u = monthly_accumulated_u / monthly_counts[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
monthly_average_v = monthly_accumulated_v / monthly_counts[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]

# Create plots for the monthly averages
for month in range(12):
    u_month = monthly_average_u[month, 0, 0, :, :]  # pressure level 0 for the u component
    v_month = monthly_average_v[month, 0, 0, :, :]  # pressure level 0 for the v component

    missing_values = np.where(u_month == nc.default_fillvals['f4'])
    print(f"Missing values in u for month {month + 1}:")
    print(u_month[missing_values])

    missing_values = np.where(v_month == nc.default_fillvals['f4'])
    print(f"Missing values in v for month {month + 1}:")
    print(v_month[missing_values])
    
    plot_spatial_map_direction(u_month, v_month, latitude, longitude, month + 1)
    # Save the figure as an image
    file_name = f"spatial_maps/monthly_map_{month+1}.png"
    plt.savefig(file_name, dpi=300)
    plt.close()

# Create a list to store the image paths
image_paths = [f"spatial_maps/monthly_map_{month+1}.png" for month in range(12)]

# Create a GIF from the images
gif_path = "spatial_maps/monthly_maps.gif"
with Image.open(image_paths[0]) as first_image:
    first_image.save(gif_path, save_all=True, append_images=[Image.open(image_path) for image_path in image_paths[1:]],
                    optimize=False, duration=200, loop=0)

# Print the size of the GIF
gif_size = os.path.getsize(gif_path)
print(f"The size of the GIF is: {gif_size} bytes")

# #not calculating monthly average version plot_spatial_map_direction_single
# def plot_spatial_map_direction_single(u, v, latitude, longitude, time):
#     wind_speed = np.sqrt(u**2 + v**2)
#     wind_direction = np.arctan2(v, u) * 180 / np.pi

#     fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
#     ax.set_extent([-162, -150, 6, 12], ccrs.PlateCarree())
#     ax.coastlines(linewidth=0.5)
#     ax.add_feature(cfeature.BORDERS, linewidth=0.5)
#     ax.add_feature(cfeature.STATES, linewidth=0.3)
#     ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

#     # Draw the region of interest box
#     lon_min, lat_min = -157, 8
#     lon_max, lat_max = -155, 10
#     ax.plot([lon_min, lon_min, lon_max, lon_max, lon_min], [lat_min, lat_max, lat_max, lat_min, lat_min],
#             transform=ccrs.PlateCarree(), color='red', linewidth=2)

#     # Add a label to the box
#     label = '*Region of\n Interest'
#     ax.text(-156, 9, label, fontsize=12, color='black', ha='center', va='center', transform=ccrs.PlateCarree())

#     # Add latitude and longitude gridlines on the bottom and left side
#     ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

#     # Plot wind speed direction with thinner arrow thickness
#     im = ax.quiver(longitude, latitude, u, v, wind_speed, cmap='jet', transform=ccrs.PlateCarree(),
#                    linewidth=0.2, clim=[np.min(wind_speed), np.max(wind_speed)])
#     # Set the arrow angles to match wind direction
#     im.set_UVC(wind_speed * np.cos(wind_direction), wind_speed * np.sin(wind_direction))

#     # Remove right and top spines
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)

#     # Add color bar (horizontal)
#     cbar = plt.colorbar(im, ax=ax, orientation='horizontal', shrink=0.7)
#     cbar.set_label('Wind Speed (m/s)', fontsize=12)

#     # Include the time in the plot title
#     plt.title(f'Wind Speed&Direction\nTime: {time}', fontsize=12)


# # Open the netCDF file
# nc_file = nc.Dataset('pacific_ocean_large.nc')
# latitude = nc_file.variables['latitude'][:]
# longitude = nc_file.variables['longitude'][:]
# u = nc_file.variables['u'][:]
# v = nc_file.variables['v'][:]
# times = nc_file.variables['time'][:]  # Get all time slots

# # Convert time to a human-readable format
# base_date = datetime.datetime(1900, 1, 1)
# converted_times = [base_date + datetime.timedelta(hours=int(time)) for time in times]

# # Create a directory to store the images
# os.makedirs('spatial_maps', exist_ok=True)
# month_indices = [i for i, t in enumerate(converted_times) if t.month <= 12]
# for i in month_indices:
#     u_time = u[i, 0, 0, :, :]
#     v_time = v[i, 0, 0, :, :]
#     converted_time = converted_times[i]

#     # Generate the spatial map for the current time slot
#     plot_spatial_map_direction_single(u_time, v_time, latitude, longitude, converted_time)

#     # Save the figure as an image
#     file_name = f"spatial_maps/spatial_map_{i}.png"
#     # file_name = "spatial_maps/spatial_map_{}.png".format(i)

#     plt.savefig(file_name, dpi=300)
#     plt.close()

# # Create a list to store the image paths
# image_paths = []
# for i in range(len(times)):
#     image_path = "spatial_maps/spatial_map_{}.png".format(i)
#     image_paths.append(image_path)

# # Create a GIF from the images
# gif_path = f"spatial_maps/spatial_map_{i}.png"
# with Image.open(image_paths[0]) as first_image:
#     first_image.save(gif_path, save_all=True, append_images=[Image.open(image_path) for image_path in image_paths[1:]],
#                     optimize=False, duration=200, loop=0)

# # Print the size of the GIF
# gif_size = os.path.getsize(gif_path)
# print(f"The size of the GIF is: {gif_size} bytes")
