import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
# Open the netCDF file
nc_file = nc.Dataset('pacific_ocean.nc')
latitude = nc_file.variables['latitude'][:]
print('latitude',latitude)
longitude = nc_file.variables['longitude'][:]
print('longitude',longitude)
u = nc_file.variables['u'][:]
v = nc_file.variables['v'][:]
wind_speed = np.sqrt(u**2 + v**2)

time = nc_file.variables['time'][:]
base_date = datetime.datetime(1900, 1, 1)
converted_dates = [base_date + datetime.timedelta(hours=int(t)) for t in time]
formatted_dates = [date.strftime('%y-%m-%d-%H-%M-%S') for date in converted_dates]
start_year = 2013
end_year = 2023

# Calculate yearly average of wind speed
yearly_average = []
for year in range(start_year, end_year + 1):
    print('year', year)  # 2013
    year_indices = [i for i, date in enumerate(converted_dates) if date.year == year]
    yearly_average.append(np.mean(wind_speed[year_indices], axis=0))
yearly_average = np.array(yearly_average)

years = np.arange(start_year, end_year + 1)
fig, ax = plt.subplots(figsize=(8, 6))
for level in range(3):
    print('level', level)
    print('nc_file.variables["level"][level]',nc_file.variables["level"][level])
    mean_yearly_average = np.mean(yearly_average[:, 0, level, :, :], axis=(1, 2))
    std_yearly_average = np.std(yearly_average[:, 0, level, :, :], axis=(1, 2))
    ax.plot(years, mean_yearly_average, '-o', label=f'Level {level} ({nc_file.variables["level"][level]} hPa)')
    ax.fill_between(years, mean_yearly_average - std_yearly_average, mean_yearly_average + std_yearly_average, alpha=0.3)

ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('Wind Speed (m/s)', fontsize=12)
ax.set_title('Yearly Average Wind Speed (2013-2023)', fontsize=14)
ax.legend()
plt.savefig('out_pic/yearly_average_ws.png')
plt.show()

hourly_average = []
for hour in range(24):
    hour_indices = [i for i, date in enumerate(converted_dates) if date.hour == hour and start_year <= date.year <= end_year]
    hourly_average.append(np.mean(wind_speed[hour_indices], axis=0))
hourly_average = np.array(hourly_average)
print('hourly_average.shape', hourly_average.shape)
# Plotting the hourly average wind speed
fig, ax = plt.subplots(figsize=(8, 6))
for level in range(3):
    print('level', level)
    print('nc_file.variables["level"][level]',nc_file.variables["level"][level])
    mean_hourly_average = np.mean(hourly_average[:,  0, level, :, :], axis=(1, 2))
    std_hourly_average = np.std(hourly_average[:,  0, level, :, :], axis=(1, 2))
    ax.plot(np.arange(24), mean_hourly_average, '-o', label=f'Level {level} ({nc_file.variables["level"][level]} hPa)')
    ax.fill_between(np.arange(24), mean_hourly_average - std_hourly_average, mean_hourly_average + std_hourly_average, alpha=0.3)

ax.set_xlabel('Hour', fontsize=14)
ax.set_ylabel('Wind Speed (m/s)', fontsize=12)
ax.set_title('Hourly Average Wind Speed (2013-2023)', fontsize=14)
ax.legend(loc='upper left', fontsize=10)
ax.tick_params(axis='x', labelsize=10)

# Format x-axis labels
hours = [datetime.time(hour=i).strftime('%H:%M:%S') for i in range(0, 24, 3)]
ax.set_xticks(range(0, 24, 3))
ax.set_xticklabels(hours, rotation=45)
plt.savefig('out_pic/hourly_average_ws.png', dpi=600)
plt.show()


#I will deal with the spatial map later!!
# Plotting the spatial map of wind speed
# plt.pcolormesh(longitude, latitude, wind_speed[0, 0, 0, :, :])  # Change the indices as per your requirement
# plt.colorbar(label='Wind Speed (m/s)')
# plt.xlabel('Longitude (degrees_east)')
# plt.ylabel('Latitude (degrees_north)')
# plt.title('Spatial Map of Wind Speed')
# plt.show()

