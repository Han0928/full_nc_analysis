import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

# Create the figure and axes
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.LambertConformal()})

# Set the extent to North America
ax.set_extent([-130, -60, 20, 55], ccrs.PlateCarree())

# Add coastlines, countries, and states
ax.coastlines(linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.3)

# Draw North America outline
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.add_feature(cfeature.LAND, facecolor='lightgray')

# Define coordinates for the box (Colorado region)
lon_min, lat_min = -109, 36
lon_max, lat_max = -102, 41

# Draw the box
ax.plot([lon_min, lon_min, lon_max, lon_max, lon_min], [lat_min, lat_max, lat_max, lat_min, lat_min],
        transform=ccrs.PlateCarree(), color='red', linewidth=2)

# Add a label to the box
label = 'Colorado'
ax.text(-104, 39, label, fontsize=8, color='black', ha='center', va='center', transform=ccrs.PlateCarree())

# Add latitude and longitude labels for the whole region
# ax.text(-125, 24, 'Lat: 20', fontsize=8, color='black', ha='right', va='bottom', transform=ccrs.PlateCarree())
# ax.text(-125, 54, 'Lat: 55', fontsize=8, color='black', ha='right', va='top', transform=ccrs.PlateCarree())
# ax.text(-58, 24, 'Lon: -130', fontsize=8, color='black', ha='left', va='bottom', transform=ccrs.PlateCarree())
# ax.text(-58, 54, 'Lon: -60', fontsize=8, color='black', ha='left', va='top', transform=ccrs.PlateCarree())

# Add topography overlay using Google Maps tiles
tiler = cimgt.GoogleTiles()
ax.add_image(tiler, 8)

# Add latitude and longitude gridlines
ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')

# Set the title
plt.title('North America \nSimulation domain: Colorado Region')
plt.savefig('output_fig/us_map.png', dpi=800)
# Display the map
plt.show()

