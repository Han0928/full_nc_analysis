import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import matplotlib.ticker as mtick
from matplotlib.dates import DateFormatter, DayLocator
import os

#test the all-day version;
# array-(list-list)-array
def read_data(file_dir):
    smps_integrated_list = []
    smps_time_datetime_list = [] 
    
    for filename in os.listdir(file_dir):
        print('filename:', filename)
        if filename.endswith('.ict'):
            print('Reading file:', filename)
            with open(os.path.join(file_dir, filename)) as f:
                header_lines = [f.readline() for i in range(7)] 
            
            year, month, day = map(int, header_lines[6].split(',')[0:3]) # the date is in the 6th line
            print('year, month, day:', year, month, day) 
            data = np.loadtxt(os.path.join(file_dir, filename), delimiter=',', skiprows=88)

            smps_integrated_list = np.concatenate((smps_integrated_list, data[:, 3]))
            print('smps_integrated_list.shape:', smps_integrated_list.shape)
            epoch = datetime.datetime(year, month, day)
            smps_time_datetime = [epoch + datetime.timedelta(seconds=t) for t in data[:, 0]] 
            smps_time_datetime_list = np.concatenate((smps_time_datetime_list, smps_time_datetime))          
            print('smps_time_datetime_list.shape', smps_time_datetime_list.shape)
    return smps_integrated_list, smps_time_datetime_list

def plot_data(smps_time_datetime_list, smps_integrated_list):
    fig, ax = plt.subplots()
    ax.scatter(smps_time_datetime_list, smps_integrated_list, s=1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Particle Concentration(5-300nm)\n(#/cm3)')
    ax.set_ylim([0, 3e5]) 
    # Format the x-axis date labels
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    plt.xticks(rotation=30, ha='right')  
    plt.show()

file_path = '/ocean/projects/atm200005p/ding0928/nc_file_full/WWW-AIR_1683138523904/'
smps_integrated_list, smps_time_datetime_list = read_data(file_path)

plot_data(smps_time_datetime_list, smps_integrated_list)
plt.savefig('output_fig/obser_BAO_tower.png', dpi=800)
