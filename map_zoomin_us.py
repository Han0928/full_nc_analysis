import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

# Create the figure and axes
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.LambertConformal()})

# Set the extent to the specified region
ax.set_extent([-110.00, -100.95, 37.95, 45.05], ccrs.PlateCarree())

# Add coastlines, countries, and states
ax.coastlines(linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.3)

# Draw North America outline
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.add_feature(cfeature.LAND, facecolor='lightgray')

# Define coordinates for the box (Colorado region)
lon_min, lat_min = -109, 36
lon_max, lat_max = -102, 42

# Draw the box as a dark gray dashed line
ax.plot([lon_min, lon_min, lon_max, lon_max, lon_min], [lat_min, lat_max, lat_max, lat_min, lat_min],
        transform=ccrs.PlateCarree(), color='darkgray', linestyle='--', linewidth=1)

# Add a label to the box
label = 'Colorado'
ax.text(-107, 42, label, fontsize=14, color='purple', ha='center', va='center', transform=ccrs.PlateCarree())

# Add latitude and longitude labels for the specified region
ax.text(-105.00, 40.03, '*BAO tower', fontsize=12, color='red', ha='left', va='bottom', transform=ccrs.PlateCarree())
ax.text(-106.6, 40.45, '*Storm Peak Lab', fontsize=12, color='red', ha='left', va='top', transform=ccrs.PlateCarree())
ax.text(-104.99, 39.73, '*Denver city', fontsize=12, color='red', ha='left', va='bottom', transform=ccrs.PlateCarree())
ax.text(-104.82, 41.14, '*Cheyenne', fontsize=12, color='red', ha='left', va='top', transform=ccrs.PlateCarree())

# Add topography overlay using Google Maps tiles
tiler = cimgt.GoogleTiles()
ax.add_image(tiler, 8)

# Add latitude and longitude gridlines
ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', xlocs=[-110, -100], ylocs=[36, 41])


# Set the title
plt.title('Simulation Domain', pad=20)  # Add a padding to adjust the title position
plt.savefig('output_fig/us_zommin.png', dpi=800)
# Display the map
plt.show()

