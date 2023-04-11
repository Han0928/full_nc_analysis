import os
from datetime import datetime, timedelta

# Loop through all files in the directory
for file_name in os.listdir('.'):
    if file_name.endswith('.nc'):  # Check if the file is a .nc file
        # Extract the time from the file name
        time_str = file_name.split('_')[-1].split('.')[0]
        
        # Convert the time string to a float representing hours since the reference time 1970
        time_hours = float(time_str)
        
        # Convert the time to a datetime object
        time_obj = datetime.fromtimestamp(time_hours*3600)
        
        # Convert the time to a human-readable format
        human_readable_time = time_obj.strftime('%Y-%m-%d %H:%M:%S')

        # Print the file name and the corresponding human-readable time
        print(f'{file_name}: {human_readable_time}')

