import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from datetime import datetime

df = pd.read_csv('/ocean/projects/atm200005p/ding0928/nc_file_full/Storm_Peak_julaug_2014.nas', 
                 delim_whitespace=True, 
                 skiprows=74, 
                 error_bad_lines=False)
ref_date = pd.to_datetime('2014-01-01')
df['date'] = ref_date + pd.to_timedelta(df['starttime'], unit='D')
df['date'] = df['date'].dt.strftime('%Y%m%d%H%M%S')
print(df['date'])

fig, ax = plt.subplots()
ax.plot(df['date'], df['particle_number_concentration'], 'bo', markersize=1.0)
ax.plot(df['date'], df['particle_number_concentration'], ':', color='blue', linewidth=0.5)
ax.set_xlabel('Date (YYYYMMDD)')
ax.set_ylabel('Particle Number Concentration \n(#/cm3)')
ax.set_title('Storm Peak Laboratory Station\n(lat: 40.445N, lon: -106.74W, alt: 3220.0 m)')

locator = mdates.DayLocator(interval=48) #2days
ax.xaxis.set_major_locator(locator)
plt.xticks(rotation=30)
ax.set_yscale('log')
ax.set_ylim([3e2, 5e4])
ax.set_xlim(pd.Timestamp('2014-07-20').strftime('%Y%m%d%H%M%S'), 
            pd.Timestamp('2014-08-10').strftime('%Y%m%d%H%M%S'))


#Now I want to add a shaded region for the top 3 days of maximum daily mean number concentration, under testing
# convert 'date' column to DatetimeIndex
df = df.set_index(pd.DatetimeIndex(df['date']))
daily_mean = daily_mean.loc['2014-07-20':'2014-08-10']
# get the top 3 days of maximum daily mean number concentration
top3_dates = daily_mean.nlargest(3).index
# convert the dates to the required format of '%Y%m%d%H%M%S'
top3_dates_str = [d.strftime('%Y%m%d') for d in top3_dates]
print(top3_dates_str)    #['20140723', '20140722', '20140803']

# add shaded regions for the top 3 days of maximum daily mean number concentration
for date_str in top3_dates_str:
    start = pd.Timestamp(date_str).strftime('%Y%m%d%H%M%S')
    end = pd.Timestamp(date_str) + pd.Timedelta(days=1)
    end = end.strftime('%Y%m%d%H%M%S')
    ax.axvspan(start, end, alpha=0.2, color='gray')
    
    # get the daily mean value for the shaded region
    daily_mean_val = daily_mean[date_str]
    daily_mean.index = pd.to_datetime(daily_mean.index.strftime('%Y%m%d'))
    
    # add an arrow showing the daily mean value
    arrow_x = pd.Timestamp(datetime.strptime(date_str, '%Y%m%d')) + pd.Timedelta(days=0.5)
    arrow_y = daily_mean_val + 1000
    ax.annotate('{:.2e}'.format(daily_mean_val), xy=(arrow_x, arrow_y), xytext=(arrow_x, arrow_y * 1.5), ha='center', arrowprops=dict(facecolor='black', arrowstyle='->'))

# the below code is the one working.
# add shaded regions for the top 3 days of maximum daily mean number concentration
# for date_str in top3_dates_str:
#     start = pd.Timestamp(date_str).strftime('%Y%m%d%H%M%S')
#     end = pd.Timestamp(date_str) + pd.Timedelta(days=1)
#     end = end.strftime('%Y%m%d%H%M%S')
#     ax.axvspan(start, end, alpha=0.2, color='gray')

plt.show()

