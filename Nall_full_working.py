# %load /ocean/projects/atm200005p/ding0928/script_full_nc/time_series_dask_working.py
import matplotlib.pyplot as plt
import os
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
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import imageio

os.environ["OPENBLAS_NUM_THREADS"] = "8"

def bbox_extract_2Dcoords(cube, bbox):
    minmax = lambda x: (np.min(x), np.max(x))
    lons = cube.coord('longitude').points 
    lats = cube.coord('latitude').points
    inregion = np.logical_and(np.logical_and(lons > bbox[0],
                                             lons < bbox[1]),
                              np.logical_and(lats > bbox[2],
                                             lats < bbox[3]))
    region_inds = np.where(inregion)
    imin, imax = minmax(region_inds[0])
    jmin, jmax = minmax(region_inds[1])
    return cube[..., imin:imax+1, jmin:jmax+1]
def add_lat_lon(cube, bbox):
    polelat = cube.coord('grid_longitude').coord_system.grid_north_pole_latitude
    polelon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude  
    source_lon = cube.coord('grid_longitude').points
    source_lat = cube.coord('grid_latitude').points
    lat2d = np.transpose(np.tile(source_lat,[len(source_lon),1]))
    lon2d = np.tile(source_lon,[len(source_lat),1])
    lons, lats = iris.analysis.cartography.unrotate_pole(lon2d, lat2d, polelon, polelat)
    longit = iris.coords.AuxCoord(lons,'longitude', units='degrees', coord_system=cs.GeogCS(6371229.0))
    latit =  iris.coords.AuxCoord(lats,'latitude', units='degrees', coord_system=cs.GeogCS(6371229.0))  
    cube.add_aux_coord(longit, (2,3)) 
    cube.add_aux_coord(latit, (2,3))
    return bbox_extract_2Dcoords(cube, bbox)
def read_pt_data(potential_temperature_file, air_pressure_file, bbox):
    potential_temperature_cube = iris.load_cube(potential_temperature_file)
    air_pressure_cube = iris.load_cube(air_pressure_file)
    potential_temperature_cube = add_lat_lon(potential_temperature_cube, bbox)
    air_pressure_cube = add_lat_lon(air_pressure_cube, bbox)
    return potential_temperature_cube, air_pressure_cube
# theta to T(k)
def convert_theta_to_temperature(potential_temperature, air_pressure):
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    Rd_cp = 287.05 / 1004.0  
    air_pressure_ratio = air_pressure/p0
    air_pressure_ratio.convert_units('1')
    temperature = potential_temperature*(air_pressure_ratio)**(Rd_cp)  
    return temperature
# from kg/kg to molecule cm-3
def mixing_ratio_to_number_concentration(mixing_ratio_data, air_pressure, actual_temperature):
    zboltz = 1.3807E-23  # (J/K) R = k * N_A, k=J/K, Avogadro's number (N_A)=6.022 x 1023 entities/mol.
    staird = air_pressure / (actual_temperature * zboltz * 1.0E6)  # 1.0E6 from m3 to cm3, another form of ideal gas law
    number_concentration = mixing_ratio_data * staird
    number_concentration.units = 'molecule cm-3'
    return number_concentration
# new version, height passed in automatically
def process_single_file(filename, air_pressure, actual_temperature, bbox, model_level=1, convert_units=True):
    if convert_units:
        variable_name = filename.split('/')[-1].split('_')[1:-1]
        variable_data_cube = iris.load_cube(filename, '_'.join(variable_name))
        variable_data_cube = add_lat_lon(variable_data_cube, bbox)
        number_concentration_data = mixing_ratio_to_number_concentration(variable_data_cube, air_pressure, actual_temperature)
        number_concentration_mean = number_concentration_data.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        number_concentration_mean = number_concentration_mean.extract(iris.Constraint(model_level_number=model_level))
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
        dia_data = add_lat_lon(dia_data, bbox)
        dia_mean = dia_data.collapsed(['grid_latitude', 'grid_longitude'], iris.analysis.MEAN)
        dia_mean = dia_mean.extract(iris.Constraint(model_level_number=model_level))
        time_data = dia_data.coord('time')
        time_data_value = time_data.points
        time_units = time_data.units
        epoch = datetime.datetime(1970, 1, 1)
        time_data_datetime = [(epoch + datetime.timedelta(hours=float(tp))) for tp in time_data_value]
        return dia_mean, time_data_datetime, model_level

def process_nc_files(filenames_N, filenames_D, air_pressure, actual_temperature, bbox):
    number_concentration_mean_values = []
    time_data_values = []
    model_level_values = []
    dia_mean_values = []
    for filename in filenames_N:
        values = process_single_file(filename, air_pressure, actual_temperature, bbox, model_level=1)
        number_concentration_mean_values.append(values[0])
        time_data_values.append(values[1])
        model_level_values.append(values[2])
        dia_mean_values.append(None) #the first 5 elements in dia_mean_values are None...
    for filename in filenames_D:
        values = process_single_file(filename, air_pressure, actual_temperature, bbox, model_level=1, convert_units=False)
        number_concentration_mean_values.append(None)
        dia_mean_values.append(values[0]) # starts from 6th element+5+5=16
        time_data_values.append(values[1])  #[1] is time_data_datetime 
        model_level_values.append(values[2]) #[2] is model_level
    return number_concentration_mean_values, time_data_values, model_level_values, dia_mean_values 

def lognormal_cumulative_forcubes(N,r,rbar,sigma):
    total=(N.data/2)*(1+sp.special.erf(np.log(r/rbar.data)/np.sqrt(2)/np.log(sigma)))
    return N.copy(total)
def cal_Nall_N100(filenames_N, filenames_D, air_pressure, actual_temperature, bbox):
    number_concentration_mean_values, time_data_values, model_level_values, dia_mean_values = process_nc_files(filenames_N, filenames_D, air_pressure, actual_temperature, bbox)
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
    print('N_100', N_100.data)
    N_all = N_nuc + N_aitken+N_acc+N_cor+N_aitin
    print('N_all', N_all.data)
    print('N<100',N_all.data-N_100.data)
    return N_all, N_100
def calculate_nucleation_contribution(cs093_data_N100, ct706_data_N100, time_data_values):
    nucleation_abso = cs093_data_N100 - ct706_data_N100
    nucleation_percentage = (nucleation_abso / ct706_data_N100) * 100
    return nucleation_abso, nucleation_percentage, time_data_values
# this is not the best way to organize the following function: one for 2 boxes
def visualize_nucleation_contribution(time_data_list, nucleation_contribution_list, bbox_list):
    try:
        fig, axes = plt.subplots(len(bbox_list), 1, figsize=(10, 8), sharex=True)
        for i, (time_data, nucleation_contribution, bbox) in enumerate(zip(time_data_list, nucleation_contribution_list, bbox_list)):
            ax = axes[i]
            print('nucleation_contribution_sorted', sorted(nucleation_contribution))
            ax.plot(time_data[i], nucleation_contribution, linestyle='--',label='Nucleation Contribution')
            ax.scatter(time_data[i], nucleation_contribution, color='blue',s=6)
            ax.axhline(0, color='black', linestyle='--')
            ax.set_xlabel('Date')
            ax.set_ylabel('(%)')
            if bbox == bbox_list[0]:
                ax.set_title('The contribution of ternary nucleation to CCN\nBAO tower')
            elif bbox == bbox_list[1]:
                ax.set_title('The contribution of ternary nucleation to CCN\nStorm Peak')
            ax.set_ylim(-60, 200)  # Set y-axis limits
            print('Domain:', bbox)
            average_value = sum(nucleation_contribution) / len(nucleation_contribution)
            print('average_value', average_value)
            ax.text(0.9, 0.8, f'Average = {average_value:.2f}%', ha='center', va='center', transform=ax.transAxes,
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))       
            ax.legend()
        plt.tight_layout()
        plt.savefig('output_fig/ccn_contri_subplot.png')
        plt.show()
    except Exception as e:
        print(f"An error occurred during nucleation contribution visualization: {str(e)}")
#one for whole domain 
def visualize_nucleation_uni(time_data_list, nucleation_contribution_list, model_level):
    fig, axes = plt.subplots(1, 1, figsize=(10, 8))  # Modify subplot creation to have 1 subplot
    ax = axes
    ax.plot(time_data_list, nucleation_contribution_list, linestyle='--', label='Nucleation Contribution')  
    ax.scatter(time_data_list, nucleation_contribution_list, color='blue',s=6) 
    ax.axhline(0, color='black', linestyle='--')
    ax.set_xlabel('Date',fontsize=14)
    ax.set_ylabel('(%)',fontsize=14)
    ax.set_title('The contribution of ternary nucleation to CCN\nColorado Region',fontsize=14)  
    ax.set_ylim(-60, 200)
    average_value = sum(nucleation_contribution_list) / len(nucleation_contribution_list)
    average_value = round(average_value, 2)
    ax.text(0.1, 0.8, f'Average = {average_value}%', ha='center', va='center', transform=ax.transAxes,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'), fontsize=14)
    ax.set_title(f"The contribution(%) of ternary nucleation to CCN\nColorado Region (Model_Level={model_level})")
    # Add shaded regions for the top 3 days of maximum contribution
    df = pd.DataFrame({'date': time_data_list, 'contribution': nucleation_contribution_list})
    df['date'] = pd.to_datetime(df['date'])
    daily_mean = df.groupby(pd.Grouper(key='date', freq='D'))['contribution'].mean()
    daily_mean = daily_mean.loc[df['date'].min():df['date'].max()]
    top3_dates = daily_mean.nlargest(3).index
    top3_dates_str = [d.strftime('%Y-%m-%d') for d in top3_dates]
    for date_str in top3_dates_str:
        start_num = mdates.date2num(pd.Timestamp(date_str))
        end_num = mdates.date2num(pd.Timestamp(date_str) + pd.Timedelta(days=1))
        ax.axvspan(start_num, end_num, alpha=0.4, color='darkgray')
        ax.text(start_num + (end_num - start_num) / 2, 0.9 * ax.get_ylim()[1],
                f"{date_str}: {round(daily_mean[date_str], 1)}",
                ha='center', va='top', fontsize=14, rotation=90)
    line_legend = mlines.Line2D([], [], color='blue', linestyle='--', label='Nucleation Contribution')
    shaded_legend = mpatches.Patch(facecolor='gray', alpha=0.4, label='Top 3 days(Average %)')
    ax.legend(handles=[line_legend, shaded_legend], loc='upper left', fontsize=12)
    plt.tight_layout()
    plt.savefig(f'output_fig/ccn_contri_wholedomain_z{model_level}.png', dpi=300)
    # List of model levels or output filenames
    output_filenames = [
        'ccn_contri_wholedomain_z1.png',
        'ccn_contri_wholedomain_z2.png',
        'ccn_contri_wholedomain_z5.png',
        'ccn_contri_wholedomain_z=10.png',
        'ccn_contri_wholedomain_z=12.png',
        'ccn_contri_wholedomain_z=15.png',
        'ccn_contri_wholedomain_z=20.png',
        'ccn_contri_wholedomain_z=25.png',
        'ccn_contri_wholedomain_z=30.png',
        'ccn_contri_wholedomain_z=40.png',
        'ccn_contri_wholedomain_z45.png'
    ]
    # Create the GIF animation from the saved images
    animation_filename = 'output_fig/ccn_contri_animation.gif'
    with imageio.get_writer(animation_filename, mode='I', duration=0.9) as writer:
        for filename in output_filenames:
            image = imageio.imread(f'output_fig/{filename}')
            writer.append_data(image)

def read_data_BAO(file_dir_BAO):
    smps_integrated_list = []
    smps_time_datetime_list = []
    for filename in os.listdir(file_dir_BAO):       
        if filename.endswith('.ict'):
            with open(os.path.join(file_dir_BAO, filename)) as f:
                header_lines = [f.readline() for i in range(7)]
            year, month, day = map(int, header_lines[6].split(',')[0:3]) # the date is in the 6th line
            data = np.loadtxt(os.path.join(file_dir_BAO, filename), delimiter=',', skiprows=88)
            smps_integrated_list = np.concatenate((smps_integrated_list, data[:, 3]))
            epoch = datetime.datetime(year, month, day)
            smps_time_datetime = [epoch+ datetime.timedelta(seconds=time_str) for time_str in data[:, 0]] 
            smps_time_datetime_list = np.concatenate((smps_time_datetime_list, smps_time_datetime))
    print('smps_time_datetime_list.value', smps_time_datetime_list)
    return smps_integrated_list, smps_time_datetime_list
# since the obser=11315, model=180, interpolate; the time format is different
def plot_Nall_comparison(time_model, smps_time, observed_data, ct706_data, cs093_data):
    time_data_unix = np.array([datetime.datetime.timestamp(dt) for dt in time_model[0][0]])
    time_model_hours = time_data_unix / 3600 
    time_smps_unix = np.array([datetime.datetime.timestamp(dt) for dt in smps_time])
    time_smps_hours = time_smps_unix / 3600 
    print('time_model_hours', time_model_hours) #[390509...391046.]
    print('time_smps_hours', time_smps_hours) #[390628.3...390987.92]
    interp_ct706 = interp1d(time_model_hours,ct706_data.data,kind='nearest', bounds_error=False,fill_value=-9999)
    ct706_data_intep = interp_ct706(np.asarray(time_smps_hours))
    print('ct706_data_interp', sorted(ct706_data_intep))
    interp_cs093 = interp1d(time_model_hours,cs093_data.data,kind='nearest', bounds_error=False,fill_value=-9999)
    cs093_data_interp = interp_cs093(np.asarray(time_smps_hours))
    print('cs093_data_intep', sorted(cs093_data_interp))

    epoch = datetime.datetime(1970, 1, 1)
    smps_time_datetime = [(epoch + datetime.timedelta(hours=float(tp))) for tp in time_smps_hours]
    fig, ax = plt.subplots() 
    ax.scatter(smps_time_datetime, observed_data, label='Observed Data: tower', s=10) #3085 nano-DMA and TSI 3025 ultra-fine CPC
    ax.scatter(smps_time_datetime, ct706_data_intep, label='Default NPF: ct706', s=10)
    ax.scatter(smps_time_datetime, cs093_data_interp, label='Updated NPF: cs093', s=10)
    for date in smps_time_datetime:
        day_start = date.replace(hour=20, minute=0, second=0)
        day_end = date.replace(hour=5, minute=30, second=0) + pd.DateOffset(days=1)
        ax.axvspan(day_start, day_end, facecolor='lightgray', alpha=0.01)
    plt.xlabel('Date', fontsize=12)
    plt.ylabel('N ($\mathrm{cm^{-3}}$)',fontsize=12)
    plt.title(' N<100nm between model and observation\n BAO tower', fontsize=12)
    plt.ylim(0, 1e5)
    plt.xticks(rotation=45)
    plt.legend(fontsize=12)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.yaxis.offsetText.set_fontsize(12)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.savefig('output_fig/Nall_com_BAO.png')

def plot_data(smps_time_datetime_list, smps_integrated_list):
    fig, ax = plt.subplots(figsize=(15,5))
    ax.scatter(smps_time_datetime_list, smps_integrated_list, s=1, color='blue', alpha=0.3)   
    ax.set_xlabel('Date', fontsize=14)
    ax.set_ylabel('Particle Concentration(5-300nm)\n(N: #/cm3)',fontsize=14)
    ax.set_ylim([0, 2e5])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))  
    ax.set_title('BAO Tower')
    ax.tick_params(axis='y', labelsize=14)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    ax.tick_params(axis='x', labelsize=14)
    plt.xticks(rotation=30, ha='right')
    print('smps_time_datetime_list',type(smps_time_datetime_list)) # <class 'numpy.ndarray'>
    # Add shaded regions for the top 3 days of maximum daily mean number concentration
    df = pd.DataFrame({'date': smps_time_datetime_list, 'particle_number_concentration': smps_integrated_list}) #pandas DataFrame df
    df['date'] = pd.to_datetime(df['date'])# datetime format 
    print('df[date]',type(df['date'])) #df[date] <class 'pandas.core.series.Series'>
    daily_mean = df.groupby(pd.Grouper(key='date', freq='D'))['particle_number_concentration'].mean() #"date" column of the df DataFrame.
    daily_mean = daily_mean.loc[df['date'].min():df['date'].max()]
    top3_dates = daily_mean.nlargest(3).index
    print('top3_dates',top3_dates)
    top3_dates_str = [d.strftime('%Y-%m-%d') for d in top3_dates] # from datetime format to string format 
    print('top3_dates_str',top3_dates_str)    # ['20140730', '20140801', '20140731']
    for date_str in top3_dates_str:
        start_num = matplotlib.dates.date2num(pd.Timestamp(date_str))
        end_num = matplotlib.dates.date2num(pd.Timestamp(date_str) + pd.Timedelta(days=1))
        ax.axvspan(start_num, end_num, alpha=0.2, color='gray')
        ax.text(start_num + (end_num - start_num) / 2, 1.3e5, f"{date_str}: N: {daily_mean[date_str]:.2e}", ha='center', va='center', fontsize=12, rotation=90)
        print('start:', start_num)
        print('end:', end_num)

def plot_data_submode(time_data_values, number_concentration_mean_values, filenames_N):
    fig, axes = plt.subplots(5, 1, figsize=(10, 20), sharex=True)
    colors = ['tab:blue', 'tab:orange']
    markers = ['o', 's']
    labels = ['Binary nucleation', 'Updated ion-ternary nucleation'] #ct076=binary; cs093=ion-ternary
    for i in range(5):
        axes[i].plot(time_data_values[i], number_concentration_mean_values[i].data, label=labels[0], color=colors[0], marker=markers[0], markersize=3, linewidth=1)
        axes[i].plot(time_data_values[i+5], number_concentration_mean_values[i+5].data, label=labels[1], color=colors[1], marker=markers[1], markersize=3, linewidth=1)
        variable_name = filenames_N[i].split('/')[-1].split('_')[8:10]
        variable_name = '_'.join(variable_name)
        title = filenames_N[i].split('/')[-1].split('_')[3].capitalize() + ' (' + variable_name + ')'
        axes[i].set_title(title, fontsize=14)
        axes[i].tick_params(axis='both', labelsize=12)
        axes[i].set_xlabel('Time', fontsize=12)
        axes[i].set_ylabel('#/cm3', fontsize=12)
        if i == 0:
            axes[i].legend()
    fig.suptitle('BAO tower,40.03(N),105.00(W) \n 0.1degree, z=1', fontsize=16, fontweight='bold') #https://psl.noaa.gov/technology/bao/site/ 
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.xticks(rotation=30)
    plt.show()
    plt.savefig('output_fig/ncfull_BAO_0.1degree_Z1.png')

# define the path of the nc files
path_ct706 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-ct706/full_nc_files/" #i_nuc=2,ion_2? i_nuc_method=2
path_cs093 = "/ocean/projects/atm200005p/ding0928/nc_file_full/u-cs093/full_nc_files/" #i_nuc=4,vn11.9_lightning_nh3npf_han

filenames_N = [
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
    path_ct706 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_insoluble_aitken_mode_aerosol_in_air_m01s34i119.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_accumulation_mode_aerosol_in_air_m01s34i107.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_aitken_mode_aerosol_in_air_m01s34i103.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_coarse_mode_aerosol_in_air_m01s34i113.nc',
    path_cs093 + 'Rgn_number_of_particles_per_air_molecule_of_soluble_nucleation_mode_aerosol_in_air_m01s34i101.nc'
]

filenames_D = [
    path_ct706 + 'Rgn_m01s38i401_m01s38i401.nc', #nuc_diam
    path_ct706 + 'Rgn_m01s38i402_m01s38i402.nc', #ait_diam
    path_ct706 + 'Rgn_m01s38i403_m01s38i403.nc', #acc_diam
    path_ct706 + 'Rgn_m01s38i404_m01s38i404.nc', #coar_diam
    path_ct706 + 'Rgn_m01s38i405_m01s38i405.nc', #ait_diam_inso
    path_cs093 + 'Rgn_m01s38i401_m01s38i401.nc',
    path_cs093 + 'Rgn_m01s38i402_m01s38i402.nc',
    path_cs093 + 'Rgn_m01s38i403_m01s38i403.nc',
    path_cs093 + 'Rgn_m01s38i404_m01s38i404.nc',
    path_cs093 + 'Rgn_m01s38i405_m01s38i405.nc'
]

potential_temperature_file_ct706 = path_ct706 + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file_ct706 = path_ct706 + 'Rgn_air_pressure_m01s00i408.nc'
potential_temperature_file_cs093 = path_cs093 + 'Rgn_air_potential_temperature_m01s00i004.nc'
air_pressure_file_cs093 = path_cs093 + 'Rgn_air_pressure_m01s00i408.nc'
box_i=1 #an index, 0 for large region, 1 for BAO tower+Storm Peak Lab
if box_i==1:
    #BAO tower,40 03(N) Longitude: 105 00(W),Storm Peak Lab:40.45° N, 106.6° W,
    bbox_list = [[-105.05, -104.95, 39.95, 40.05], [-106.65, -106.55, 40.40, 40.50]]  # dirty first, clean later
    time_data_list = [] 
    nucleation_percentage_list = [] 
    # bbox_list[0]
    potential_temperature_ct706, air_pressure_ct706 = read_pt_data(potential_temperature_file_ct706, air_pressure_file_ct706, bbox_list[0])
    actual_temperature_ct706 = convert_theta_to_temperature(potential_temperature_ct706, air_pressure_ct706)
    potential_temperature_cs093, air_pressure_cs093 = read_pt_data(potential_temperature_file_cs093, air_pressure_file_cs093, bbox_list[0])  # Use bbox_list[0] for both CT706 and CS093
    actual_temperature_cs093 = convert_theta_to_temperature(potential_temperature_cs093, air_pressure_cs093)

    number_concentration_mean_values_ct706, time_data_values_ct706, model_level_values_ct706, dia_mean_values_ct706 = process_nc_files(filenames_N[:5], filenames_D[:5], air_pressure_ct706, actual_temperature_ct706, bbox_list[0])
    number_concentration_mean_values_cs093, time_data_values_cs093, model_level_values_cs093, dia_mean_values_cs093 = process_nc_files(filenames_N[5:], filenames_D[5:], air_pressure_cs093, actual_temperature_cs093, bbox_list[0])

    n_all_mixing_ct706, n100_mixing_ct706 = cal_Nall_N100(filenames_N[:5], filenames_D[:5], air_pressure_ct706, actual_temperature_ct706, bbox_list[0])
    n_all_mixing_cs093, n100_mixing_cs093 = cal_Nall_N100(filenames_N[5:], filenames_D[5:], air_pressure_cs093, actual_temperature_cs093, bbox_list[0])
    nucleation_abso, nucleation_percentage, time_data = calculate_nucleation_contribution(n100_mixing_cs093, n100_mixing_ct706, time_data_values_ct706)
    time_data_list.append(time_data)
    nucleation_percentage_list.append(nucleation_percentage.data)

    file_path_BAO = '/ocean/projects/atm200005p/ding0928/nc_file_full/WWW-AIR_1683138523904/'
    smps_integrated_list, smps_time_datetime_list = read_data_BAO(file_path_BAO)
    plot_Nall_comparison(time_data_list, smps_time_datetime_list, smps_integrated_list, n_all_mixing_ct706-n100_mixing_ct706, n_all_mixing_cs093-n100_mixing_cs093)

    # bbox_list[1]
    potential_temperature_ct706, air_pressure_ct706 = read_pt_data(potential_temperature_file_ct706, air_pressure_file_ct706, bbox_list[1])
    actual_temperature_ct706 = convert_theta_to_temperature(potential_temperature_ct706, air_pressure_ct706)
    potential_temperature_cs093, air_pressure_cs093 = read_pt_data(potential_temperature_file_cs093, air_pressure_file_cs093, bbox_list[1])
    actual_temperature_cs093 = convert_theta_to_temperature(potential_temperature_cs093, air_pressure_cs093)

    number_concentration_mean_values_ct706, time_data_values_ct706, model_level_values_ct706, dia_mean_values_ct706 = process_nc_files(filenames_N[:5], filenames_D[:5], air_pressure_ct706, actual_temperature_ct706, bbox_list[1])
    number_concentration_mean_values_cs093, time_data_values_cs093, model_level_values_cs093, dia_mean_values_cs093 = process_nc_files(filenames_N[5:], filenames_D[5:], air_pressure_cs093, actual_temperature_cs093, bbox_list[1])

    n_all_mixing_ct706, n100_mixing_ct706 = cal_Nall_N100(filenames_N[:5], filenames_D[:5], air_pressure_ct706, actual_temperature_ct706, bbox_list[1])
    n_all_mixing_cs093, n100_mixing_cs093 = cal_Nall_N100(filenames_N[5:], filenames_D[5:], air_pressure_cs093, actual_temperature_cs093, bbox_list[1])
    nucleation_abso, nucleation_percentage, time_data = calculate_nucleation_contribution(n100_mixing_cs093, n100_mixing_ct706, time_data_values_ct706) # time_data_values_ct706 or cs093, doesn't matter

    time_data_list.append(time_data)
    nucleation_percentage_list.append(nucleation_percentage.data)
    visualize_nucleation_contribution(time_data_list, nucleation_percentage_list, bbox_list)  
else:
    bbox = [-108.00, -102.95, 37.95, 43.05] #whole domain,(40.0,-105)
    potential_temperature_ct706, air_pressure_ct706 = read_pt_data(potential_temperature_file_ct706, air_pressure_file_ct706, bbox)
    actual_temperature_ct706 = convert_theta_to_temperature(potential_temperature_ct706, air_pressure_ct706)

    potential_temperature_cs093, air_pressure_cs093 = read_pt_data(potential_temperature_file_cs093, air_pressure_file_cs093, bbox)
    actual_temperature_cs093 = convert_theta_to_temperature(potential_temperature_cs093, air_pressure_cs093)

    number_concentration_mean_values_ct706, time_data_values_ct706,model_level_values_ct706, dia_mean_values_ct706 = process_nc_files(filenames_N[:5], filenames_D[:5],air_pressure_ct706, actual_temperature_ct706, bbox)
    number_concentration_mean_values_cs093, time_data_values_cs093,model_level_values_cs093, dia_mean_values_cs093 = process_nc_files(filenames_N[5:], filenames_D[5:],air_pressure_cs093, actual_temperature_cs093, bbox)

    n_all_mixing_ct706, n100_mixing_ct706 = cal_Nall_N100(filenames_N[:5], filenames_D[:5], air_pressure_ct706, actual_temperature_ct706, bbox)
    n_all_mixing_cs093, n100_mixing_cs093 = cal_Nall_N100(filenames_N[5:], filenames_D[5:], air_pressure_ct706, actual_temperature_ct706, bbox)
    nucleation_abso, nucleation_percentage, time_data = calculate_nucleation_contribution(n100_mixing_cs093, n100_mixing_ct706, time_data_values_ct706)
    print('model_level_values_cs093[0]_value',model_level_values_cs093[0])
    visualize_nucleation_uni(time_data[0], nucleation_percentage.data, model_level=model_level_values_cs093[0])

    number_concentration_mean_values = number_concentration_mean_values_ct706 + number_concentration_mean_values_cs093
    time_data_values = time_data_values_ct706 + time_data_values_cs093
    
    # plot_data_submode(time_data_values, number_concentration_mean_values, filenames_N[:5] + filenames_N[5:])
    # plot_data_submode(time_data, number_concentration_mean_values, filenames_N[:5] + filenames_N[5:])


#some understanding for myself of the loop:
# so,the loop order, for ct706, the ending list is 10 for each: 
# 247-79/93-127/128, for filenames_N, executeed 5 times;numb_list[# cc-3,# cc-3,# cc-3,# cc-3,# cc-3,NA,NA,NA,NA,NA]
# 247-97-133/134, for filenames_D;[NA,NA,NA,NA,NA,# cc-3,# cc-3,# cc-3,# cc-3,# cc-3]

# then, for cs093, repeat the same loop again?
# 248-79/93-127/128, for filenames_N, executeed 5 times;numb_list[# cc-3,# cc-3,# cc-3,# cc-3,# cc-3,NA,NA,NA,NA,NA]
# 248-97-133/134, for filenames_D;[NA,NA,NA,NA,NA,# cc-3,# cc-3,# cc-3,# cc-3,# cc-3]
# the above loop makes sense now.
# then cal_Nall_N100,
# 253-79/93-128/129, for filenames_N, executeed 5 times;numb_list[# cc-3,# cc-3,# cc-3,# cc-3,# cc-3,NA,NA,NA,NA,NA]
# 253-97-135/136, for filenames_D;[NA,NA,NA,NA,NA,# cc-3,# cc-3,# cc-3,# cc-3,# cc-3]
# then print out,15 lines