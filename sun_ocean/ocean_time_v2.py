import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from PIL import Image

def calculate_yearly_average(wind_speed, converted_dates, start_year, end_year):
    yearly_average = []
    for year in range(start_year, end_year + 1):
        year_indices = [i for i, date in enumerate(converted_dates) if date.year == year]
        yearly_average.append(np.mean(wind_speed[year_indices], axis=0))
    return np.array(yearly_average)

def calculate_hourly_average(wind_speed, converted_dates, start_year, end_year):
    hourly_average = []
    for hour in range(24):
        hour_indices = [i for i, date in enumerate(converted_dates) if date.hour == hour and start_year <= date.year <= end_year]
        hourly_average.append(np.mean(wind_speed[hour_indices], axis=0))
    return np.array(hourly_average)

def calculate_monthly_mean(wind_speed, converted_dates, start_year, end_year):
    monthly_mean = []
    for month in range(1, 13):
        month_indices = [
            i for i, date in enumerate(converted_dates)
            if date.month == month and start_year <= date.year <= end_year
        ]
        monthly_mean.append(np.mean(wind_speed[month_indices], axis=0))
    return np.array(monthly_mean)

def plot_yearly_average(years, yearly_average, levels, level_names):
    fig, ax = plt.subplots(figsize=(8, 6))
    for level in range(levels):
        mean_yearly_average = np.mean(yearly_average[:, 0, level, :, :], axis=(1, 2))
        std_yearly_average = np.std(yearly_average[:, 0, level, :, :], axis=(1, 2))
        ax.plot(years, mean_yearly_average, '-o', label=f'Level {level} ({level_names[level]} hPa)')
        ax.fill_between(years, mean_yearly_average - std_yearly_average, mean_yearly_average + std_yearly_average, alpha=0.3)

    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Wind Speed (m/s)', fontsize=12)
    ax.set_title('Yearly Average Wind Speed (2013-2023)', fontsize=14)
    ax.legend()
    plt.savefig('out_pic/yearly_average_ws.png', dpi=600)
    # plt.show()

def plot_hourly_average(hourly_average, levels, level_names):
    fig, ax = plt.subplots(figsize=(8, 6))
    for level in range(levels):
        mean_hourly_average = np.mean(hourly_average[:,  0, level, :, :], axis=(1, 2))
        std_hourly_average = np.std(hourly_average[:,  0, level, :, :], axis=(1, 2))
        ax.plot(np.arange(24), mean_hourly_average, '-o', label=f'Level {level} ({level_names[level]} hPa)')
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

def plot_monthly_mean(monthly_mean, levels, level_names):
    fig, ax = plt.subplots(figsize=(8, 6))
    for level in range(levels):
        mean_monthly_mean = np.mean(monthly_mean[:, 0, level, :, :], axis=(1, 2))
        std_monthly_mean = np.std(monthly_mean[:, 0, level, :, :], axis=(1, 2))
        ax.plot(range(1, 13), mean_monthly_mean, '-o', label=f'Level {level} ({level_names[level]} hPa)')
        ax.fill_between(range(1, 13), mean_monthly_mean - std_monthly_mean, mean_monthly_mean + std_monthly_mean, alpha=0.3)

    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Wind Speed (m/s)', fontsize=12)
    ax.set_title('Monthly Mean Wind Speed (2013-2023)', fontsize=14)
    ax.legend()
    # Format x-axis labels
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(months)
    plt.savefig('out_pic/monthly_mean_ws.png', dpi=600)

def process_wind_speed_data(nc_file_path, start_year, end_year):
    # Open the netCDF file
    nc_file = nc.Dataset(nc_file_path)
    latitude = nc_file.variables['latitude'][:]
    longitude = nc_file.variables['longitude'][:]
    u = nc_file.variables['u'][:]
    v = nc_file.variables['v'][:]
    # wind_speed = np.sqrt(u**2 + v**2)
    wind_speed = np.sqrt(np.square(u) + np.square(v))
    wind_speed = np.ma.masked_invalid(wind_speed)

    time = nc_file.variables['time'][:]
    base_date = datetime.datetime(1900, 1, 1)
    converted_dates = [base_date + datetime.timedelta(hours=int(t)) for t in time]

    # Calculate yearly average of wind speed
    yearly_average = calculate_yearly_average(wind_speed, converted_dates, start_year, end_year)
    years = np.arange(start_year, end_year + 1)
    levels = len(nc_file.variables["level"])
    level_names = nc_file.variables["level"][:]
    plot_yearly_average(years, yearly_average, levels, level_names)

    # Calculate hourly average of wind speed
    hourly_average = calculate_hourly_average(wind_speed, converted_dates, start_year, end_year)
    levels = len(nc_file.variables["level"])
    level_names = nc_file.variables["level"][:]
    plot_hourly_average(hourly_average, levels, level_names)

    # Calculate monthly mean of wind speed
    monthly_mean = calculate_monthly_mean(wind_speed, converted_dates, start_year, end_year)
    levels = len(nc_file.variables["level"])
    level_names = nc_file.variables["level"][:]
    # Plot monthly mean wind speed
    plot_monthly_mean(monthly_mean, levels, level_names)

process_wind_speed_data('pacific_ocean.nc', 2013, 2023)
