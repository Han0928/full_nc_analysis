# %load /ocean/projects/atm200005p/ding0928/script_full_nc/time_series_dask_working.py
"""
this script is mainly dealing with the concanated.nc date located in /ocean/projects/atm200005p/ding0928/nc_file_full/ucy628/full_nc_files/
where ucy628 is the bhn scheme, and ucy627 is the ion-ternary scheme.
step 1: read in the data files, and extract the data from the files
step 2: calculate the Nall and N100
step 3: post-process the data, and visualize the data: time series
Necessary: interact -n 128 in terminal, otherwise it's too slow to run the script.
The github version is sit in: https://github.com/Han0928/full_nc_analysis.git 
-Han, 2023 August
"""
from shapely.geometry import box
import matplotlib.pyplot as plt
import iris.plot as iplt
import os
import imageio
import numpy as np
import iris
import iris.coord_systems as cs
import iris.coord_systems as coord_systems
import datetime
from scipy.interpolate import interp1d, RegularGridInterpolator
import scipy as sp
import netCDF4 as nc
import pandas as pd
import matplotlib.dates as mdates
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import imageio
from pathlib import Path
from matplotlib.patches import FancyArrowPatch
from datetime import timedelta
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature

os.environ["OPENBLAS_NUM_THREADS"] = "8"
def bbox_extract(cube, bbox):
    if bbox is None:
        raise ValueError("Bounding box (bbox) cannot be None.")
    print("Before_Extracting data for bbox:", bbox) 
    lat_constraint = iris.Constraint(latitude=lambda x: bbox[2] <= x <= bbox[3])
   
    cube = cube.extract(lat_constraint)# & lon_constraint)
    cube = cube.intersection(longitude=(bbox[0],bbox[1]))
    print("After_Extracting data for bbox:", bbox)
    if cube is None or cube.data.size == 0:
        print("No data extracted for bbox:", bbox)
    else:
        print("Successfully extracted data for bbox:", bbox)
    return cube

def read_pt_data(potential_temperature_file, air_pressure_file, bbox):
    potential_temperature_cube = iris.load_cube(potential_temperature_file)
    air_pressure_cube = iris.load_cube(air_pressure_file)
    potential_temperature_cube = bbox_extract(potential_temperature_cube, bbox)
    air_pressure_cube = bbox_extract(air_pressure_cube, bbox)
    return potential_temperature_cube, air_pressure_cube

def convert_theta_to_temperature(potential_temperature, air_pressure):
    original_name = potential_temperature.name()
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    Rd_cp = 287.05 / 1004.0  
    air_pressure_ratio = air_pressure/p0
    air_pressure_ratio.convert_units('1')
    for cube in [potential_temperature, air_pressure_ratio]:
        for factory in cube.aux_factories:
            if factory.name() == 'altitude':
                cube.remove_aux_factory(factory)
    temperature = potential_temperature*(air_pressure_ratio)**(Rd_cp)
    return temperature, original_name

def remove_altitude_coord(cube, cube_name):
    try:
        altitude_coord = cube.coord('altitude')
        cube.remove_coord(altitude_coord)
    except iris.exceptions.CoordinateNotFoundError:
        print(f"Altitude coordinate not found in {cube_name} cube.")

def mixing_ratio_to_number_concentration(mixing_ratio_data, air_pressure, actual_temperature):
    zboltz = 1.3807E-23  # (J/K) R = k * N_A, k=J/K, Avogadro's number (N_A)=6.022 x 1023 entities/mol.
    original_name = mixing_ratio_data.name()
    remove_altitude_coord(air_pressure, 'air_pressure')
    staird = air_pressure / (actual_temperature * zboltz * 1.0E6)  # 1.0E6 from m3 to cm3, another form of ideal gas law
    remove_altitude_coord(mixing_ratio_data, 'mixing_ratio_data')
    number_concentration = mixing_ratio_data * staird
    number_concentration.units = 'molecule cm-3'
    return number_concentration, original_name

def process_single_file(filename, air_pressure, actual_temperature, bbox, model_level=1, convert_units=True):
    if convert_units:
        variable_name = filename.split('/')[-1].split('_')[1:-1]
        variable_data_cube = iris.load_cube(filename, '_'.join(variable_name))
        variable_data_cube = bbox_extract(variable_data_cube, bbox)
        number_concentration_data, original_name = mixing_ratio_to_number_concentration(variable_data_cube, air_pressure, actual_temperature)
        number_concentration_data.rename(original_name) # Reapply the original name
        number_concentration_mean = number_concentration_data.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
        number_concentration_mean = number_concentration_mean.extract(iris.Constraint(atmosphere_hybrid_height_coordinate=model_level))
        time_data = variable_data_cube.coord('time')
        time_data_value = time_data.points
        time_units = time_data.units
        epoch = datetime.datetime(1970, 1, 1)
        time_data_datetime = [(epoch + datetime.timedelta(hours=float(tp))) for tp in time_data_value]  #<class 'list'>
        return number_concentration_mean, time_data_datetime, model_level
    else:
        with nc.Dataset(filename, 'r') as ncfile:
            variable_name = list(ncfile.variables.keys())[0]
        dia_data = iris.load_cube(filename, variable_name)
        dia_data = bbox_extract(dia_data, bbox)
        dia_mean = dia_data.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
        dia_mean = dia_mean.extract(iris.Constraint(model_level_number=model_level))
        time_data = dia_data.coord('time')
        time_data_value = time_data.points
        time_units = time_data.units
        epoch = datetime.datetime(1970, 1, 1)
        time_data_datetime = [(epoch + datetime.timedelta(hours=float(tp))) for tp in time_data_value]
        print("time_data_datetime: ", time_data_datetime)
        return dia_mean, time_data_datetime, model_level

def process_nc_files(filenames_N, filenames_D, air_pressure, actual_temperature, bbox, model_level):
    number_concentration_mean_values = []
    time_data_values = []
    model_level_values = []
    dia_mean_values = []
    for filename in filenames_N:
        values = process_single_file(filename, air_pressure, actual_temperature, bbox, model_level=model_level)
        number_concentration_mean_values.append(values[0])
        time_data_values.append(values[1])
        model_level_values.append(values[2])
        dia_mean_values.append(None) #the first 5 elements in dia_mean_values are None...
    for filename in filenames_D:
        values = process_single_file(filename, air_pressure, actual_temperature, bbox, model_level=model_level, convert_units=False)
        number_concentration_mean_values.append(None)
        dia_mean_values.append(values[0]) # starts from 6th element+5+5=16
        time_data_values.append(values[1])  #[1] is time_data_datetime 
        model_level_values.append(values[2]) #[2] is model_level
    return number_concentration_mean_values, time_data_values, model_level_values, dia_mean_values 

# The accumulation function from the lower end to the upper end of the distribution;
def lognormal_cumulative_forcubes(N,r,rbar,sigma):
    total=(N.data/2)*(1+sp.special.erf(np.log(r/rbar.data)/np.sqrt(2)/np.log(sigma)))
    return N.copy(total)

def calculate_Nall_N100(filenames_N, filenames_D, air_pressure, actual_temperature, bbox):
    number_concentration_mean_values, time_data_values, model_level_values, dia_mean_values = process_nc_files(filenames_N, filenames_D, air_pressure, actual_temperature, bbox, model_level)
    # pay very much attention to the order of the following 5 lines, depends on the order of you read in the N files #210.
    N_nuc, D_nuc = number_concentration_mean_values[4], dia_mean_values[5] # nucleation mode
    N_aitken, D_aitken = number_concentration_mean_values[2], dia_mean_values[6] # Aitken mode
    N_acc, D_acc = number_concentration_mean_values[1], dia_mean_values[7] # accumulation mode
    N_cor, D_cor = number_concentration_mean_values[3], dia_mean_values[8] # coarse mode
    N_aitin, D_aitin = number_concentration_mean_values[0], dia_mean_values[9] # insoluble Aitken mode
    # Calculate the nucleation/CCN contribution
    N_100= N_aitken - lognormal_cumulative_forcubes(N_aitken, 1.0e-7, D_aitken, 1.59) #1.59: soluble nucleation and Aitken modes,
    N_100 += N_aitin - lognormal_cumulative_forcubes(N_aitin, 1.0e-7, D_aitin, 1.59) #1.59:insoluble Aitken and accumulation modes,
    N_100+= N_acc - lognormal_cumulative_forcubes(N_acc, 1.0e-6, D_acc, 1.4) ##1.4:accumulation mode, # 2.0: coarse mode.
    N_100 += N_cor
    N_all = N_nuc + N_aitken+N_acc+N_cor+N_aitin
    return N_all, N_100, N_nuc

def calculate_nucleation_contribution(cy627_data_N100, cy628_data_N100, time_data_values):
    nucleation_abso = cy627_data_N100 - cy628_data_N100
    nucleation_percentage = (nucleation_abso / cy628_data_N100) * 100
    return nucleation_abso, nucleation_percentage, time_data_values

def plot_time_series(time_data_list, nucleation_contribution_list, region_names, bboxes, model_level):
    fig = plt.figure(figsize=(18, 8), dpi=180)
    ax1 = fig.add_subplot(1, 2, 1)
    for time_data, nucleation_percentage, region_name in zip(time_data_list, nucleation_contribution_list, region_names):
        ax1.plot(time_data[0], nucleation_percentage)
        # Finding the peak point on each line to place my legend(city name)
        peak_index = np.argmax(nucleation_percentage)
        peak_point_x = time_data[0][peak_index]
        peak_point_y = nucleation_percentage[peak_index]
        # Determining the slope at the peak to decide the annotation position
        if peak_index > 0 and peak_index < len(time_data[0]) - 1:
            slope = nucleation_percentage[peak_index + 1] - nucleation_percentage[peak_index - 1]
        else:
            slope = 0
        offset = (-5, 5) if slope < 0 else (5, 5)
        offset_days = offset[0]
        offset_time = timedelta(days=offset_days) # Creating a timedelta object
        new_point_x = peak_point_x + offset_time

        arrow = FancyArrowPatch((peak_point_x, peak_point_y), 
                                (new_point_x, peak_point_y + offset[1]),
                                arrowstyle='-|>', mutation_scale=15)
        plt.gca().add_patch(arrow)

        plt.text(new_point_x, peak_point_y + offset[1], region_name,
                 fontsize=12, family='serif')
        ax1.set_xlabel('Date', fontsize=14, family='serif')
        ax1.set_ylabel('Nucleation Percentage(%)', fontsize=14, family='serif')
        ax1.set_title(f'(a) Global Updated Ternary Nucleation parameterization (z={model_level})', fontsize=16, family='serif')
        ax1.tick_params(labelsize=12)
        ax1.grid(True, linestyle='--', alpha=0.5)  
    
    # Create the right subplot for the map
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.BORDERS, linestyle=':')
    
    for bbox, region_name in zip(bboxes, region_names):
        x = [bbox[0], bbox[0], bbox[1], bbox[1], bbox[0]]
        y = [bbox[2], bbox[3], bbox[3], bbox[2], bbox[2]]
        ax2.plot(x, y, 'r', transform=ccrs.PlateCarree())
        center_x = (bbox[0] + bbox[1]) / 2
        center_y = (bbox[2] + bbox[3]) / 2
        ax2.text(center_x, center_y, region_name, ha='center', fontsize=10, transform=ccrs.PlateCarree(), color='blue')
    ax2.set_title(f'(b) Region Names (z={model_level})', fontsize=16, family='serif')
    plt.tight_layout()
    plt.savefig(f"output_fig/time_z={model_level}.png", dpi=300, bbox_inches='tight')

def write_to_excel(time_data, N_all_values, N_100_values, region_names, filename="output.xlsx"):
    print("N_all_values shape:", np.shape(N_all_values))  # Add this line
    print("N_100_values shape:", np.shape(N_100_values))  # Add this line
    rows = []
    for month in range(12):
        row = {'Date': time_data[month]}
        for i, loc in enumerate(region_names):
            row[f'{loc}_N_all'] = N_all_values[i][month]
            row[f'{loc}_N_100'] = N_100_values[i][month]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_excel(filename, index=False)

# the main function 
rose_id=['gl_cy628','gl_cy627']
path_cy628 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-cy628/full_nc_files/" #i_nuc=2,ion_2? i_nuc_method=2
path_cy627 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-cy627/full_nc_files/" #i_nuc=4,vn11.9_lightning_nh3npf_han

filenames_N = [
    path_cy628 + 'ucy628_m01s34i119.nc',
    path_cy628 + 'ucy628_m01s34i107.nc',
    path_cy628 + 'ucy628_m01s34i103.nc',
    path_cy628 + 'ucy628_m01s34i113.nc',
    path_cy628 + 'ucy628_m01s34i101.nc',
    path_cy627 + 'ucy627_m01s34i119.nc',
    path_cy627 + 'ucy627_m01s34i107.nc',
    path_cy627 + 'ucy627_m01s34i103.nc',
    path_cy627 + 'ucy627_m01s34i113.nc',
    path_cy627 + 'ucy627_m01s34i101.nc'
]

filenames_D = [
    path_cy628 + 'ucy628_m01s38i401.nc', #nuc_diam
    path_cy628 + 'ucy628_m01s38i402.nc', #ait_diam
    path_cy628 + 'ucy628_m01s38i403.nc', #acc_diam
    path_cy628 + 'ucy628_m01s38i404.nc', #coar_diam
    path_cy628 + 'ucy628_m01s38i405.nc', #ait_diam_inso
    path_cy627 + 'ucy627_m01s38i401.nc',
    path_cy627 + 'ucy627_m01s38i402.nc',
    path_cy627 + 'ucy627_m01s38i403.nc',
    path_cy627 + 'ucy627_m01s38i404.nc',
    path_cy627 + 'ucy627_m01s38i405.nc'
]

potential_temperature_file_cy628 = path_cy628 + 'ucy628_m01s00i004.nc'
air_pressure_file_cy628 = path_cy628 + 'ucy628_m01s00i408.nc'
potential_temperature_file_cy627 = path_cy627 + 'ucy627_m01s00i004.nc'
air_pressure_file_cy627 = path_cy627 + 'ucy627_m01s00i408.nc'
box_i=1
if box_i==1:
    bbox_list = [
        [-140, -130, 50, 60], [-170, -160, 75, 85], [-100, -90, 40, 50],
        [-60, -50, -30, -20], [10, 20, 10, 20], [0, 10, 50, 60],
        [30, 40, 50, 60], [40, 50, 30, 40], [70, 80, 20, 30], [110, 120, 30, 40]
    ]
    region_names = [
    "Pacific Northwest (NA)", "Arctic Region",
    "Midwest USA", "South America",
    "Central Africa", "North Europe",
    "East Europe", "Central Asia",
    "South Asia", "East Asia"]

    bboxes = [[-140, -130, 50, 60], [-170, -160, 75, 85], [-100, -90, 40, 50],
        [-60, -50, -30, -20], [10, 20, 10, 20], [0, 10, 50, 60],
        [30, 40, 50, 60], [40, 50, 30, 40], [70, 80, 20, 30], [110, 120, 30, 40]]
    
    model_levels = [1, 5, 10, 15, 20, 25, 30, 40]
    png_files = []

    for model_level in model_levels:
        time_data_list = []
        nucleation_percentage_list = []
        n_nuc_mixing_cy628_list = [] # default i_nuc=2
        n_nuc_mixing_cy627_list = [] # updated i_nuc=4
        n_all_mixing_cy628_list = [] # I need this list to write to excel
        n100_mixing_cy628_list = []
        n_all_mixing_cy627_list = []
        n100_mixing_cy627_list = []
    
        for bbox in bbox_list:
            potential_temperature_cy628, air_pressure_cy628 = read_pt_data(potential_temperature_file_cy628, air_pressure_file_cy628, bbox)
            actual_temperature_cy628, original_name = convert_theta_to_temperature(potential_temperature_cy628, air_pressure_cy628)
            actual_temperature_cy628.rename(original_name)

            potential_temperature_cy627, air_pressure_cy627 = read_pt_data(potential_temperature_file_cy627, air_pressure_file_cy627, bbox)
            actual_temperature_cy627, original_name = convert_theta_to_temperature(potential_temperature_cy627, air_pressure_cy627)
            actual_temperature_cy627.rename(original_name)
            
            number_concentration_mean_values_cy628, time_data_values_cy628, model_level_values_cy628, dia_mean_values_cy628 = process_nc_files(filenames_N[:5], filenames_D[:5], air_pressure_cy628, actual_temperature_cy628, bbox, model_level)
            number_concentration_mean_values_cy627, time_data_values_cy627, model_level_values_cy627, dia_mean_values_cy627 = process_nc_files(filenames_N[5:], filenames_D[5:], air_pressure_cy627, actual_temperature_cy627, bbox, model_level)
        
            n_all_mixing_cy628, n100_mixing_cy628, n_nuc_mixing_cy628 = calculate_Nall_N100(filenames_N[:5], filenames_D[:5], air_pressure_cy628, actual_temperature_cy628, bbox)
            n_all_mixing_cy627, n100_mixing_cy627, n_nuc_mixing_cy627 = calculate_Nall_N100(filenames_N[5:], filenames_D[5:], air_pressure_cy627, actual_temperature_cy627, bbox)

            nucleation_abso, nucleation_percentage, time_data = calculate_nucleation_contribution(n100_mixing_cy627, n100_mixing_cy628, time_data_values_cy628)

            time_data_list.append(time_data)
            nucleation_percentage_list.append(nucleation_percentage.data)

            n_nuc_mixing_cy628_list.append(n_nuc_mixing_cy628.data)
            n_nuc_mixing_cy627_list.append(n_nuc_mixing_cy627.data)

            n_all_mixing_cy628_list.append(n_all_mixing_cy628.data)
            n100_mixing_cy628_list.append(n100_mixing_cy628.data)
            n_all_mixing_cy627_list.append(n_all_mixing_cy627.data)
            n100_mixing_cy627_list.append(n100_mixing_cy627.data)

        plot_time_series(time_data_list, nucleation_percentage_list, region_names, bboxes, model_level)

        png_filename = f'output_fig/time_z={model_level}.png'
        png_files.append(png_filename)
        # Excel outputs for each model level after processing all bboxes
        # write_to_excel(time_data_list[0][0], n_all_mixing_cy628_list, n100_mixing_cy628_list, region_names, filename=f"628output_level_{model_level}.xlsx")
        # write_to_excel(time_data_list[0][0], n_all_mixing_cy627_list, n100_mixing_cy627_list, region_names, filename=f"627output_level_{model_level}.xlsx")
        print(f"Finished processing for bbox {bbox} at model level {model_level}.")
    
    # Convert saved PNGs to GIF
    output_gif = "time_animation.gif"
    with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
        for png_file in png_files:
            image = imageio.imread(png_file)
            writer.append_data(image)
    print(f"GIF saved as {output_gif}")   
plt.show()

''' 
def bbox_extract(cube, bbox):
    bbox[0] = bbox[0] % 360
    bbox[1] = bbox[1] % 360
    lat_constraint = iris.Constraint(latitude=lambda x: bbox[2] <= x <= bbox[3])
    lon_constraint = iris.Constraint(longitude=lambda x: bbox[0] <= x <= bbox[1])
    cube = cube.extract(lat_constraint & lon_constraint)

    # Getting the latitude and longitude coordinates
    latitudes = cube.coord('latitude').points
    longitudes = cube.coord('longitude').points
    print("_after_Latitudes:", latitudes)
    print("_after_Longitudes:", longitudes)
    return cube
'''