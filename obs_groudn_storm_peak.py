# function version
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from datetime import datetime
import matplotlib.offsetbox as offsetbox

def read_data(file_path):
    df = pd.read_csv(file_path, delim_whitespace=True, skiprows=74, error_bad_lines=False)
    ref_date = pd.to_datetime('2014-01-01')
    df['date'] = ref_date + pd.to_timedelta(df['starttime'], unit='D')
    df['date'] = df['date'].dt.strftime('%Y%m%d%H%M%S')
    return df

def plot_data(df):
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.plot(df['date'], df['particle_number_concentration'], 'bo', markersize=1.5)
    ax.plot(df['date'], df['particle_number_concentration'], ':', color='blue', linewidth=1.5)
    ax.set_xlabel('Date (YYYYMMDD)', fontsize=14)
    ax.set_ylabel('Particle Number Concentration \n(#/cm3)', fontsize=14)
    ax.set_title('Storm Peak Laboratory Station\n(lat: 40.445N, lon: -106.74W, alt: 3220.0 m)')

    locator = mdates.DayLocator(interval=48)
    ax.xaxis.set_major_locator(locator)
    plt.xticks(rotation=30)
    ax.set_yscale('log')
    ax.set_ylim([3e2, 5e4])
    ax.set_xlim(pd.Timestamp('2014-07-20').strftime('%Y%m%d%H%M%S'), pd.Timestamp('2014-08-10').strftime('%Y%m%d%H%M%S'))
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.get_major_formatter().set_powerlimits((0, 4))
    ax.tick_params(axis='y', labelsize=12)
    return ax

def add_shaded_regions(ax, df):
    #Now I want to add a shaded region for the top 3 days of maximum daily mean number concentration
    # step 1:calculate the daily mean number concentration
    df['date'] = pd.to_datetime(df['date'])
    daily_mean = df.groupby(pd.Grouper(key='date', freq='D'))['particle_number_concentration'].mean()
    # step 2: get the top 3 days of maximum daily mean number concentration
    df = df.set_index(pd.DatetimeIndex(df['date']))
    daily_mean = daily_mean.loc['2014-07-20':'2014-08-10']
    top3_dates = daily_mean.nlargest(3).index
    top3_dates_str = [d.strftime('%Y%m%d') for d in top3_dates]
    print('top3_dates_str',top3_dates_str)    # ['20140723', '20140722', '20140803']
    # step 3: add shaded regions for the top 3 days of maximum daily mean number concentration failed!

    #version 3 for arrow:
    for date_str in top3_dates_str:
        start = pd.Timestamp(date_str).strftime('%Y%m%d%H%M%S')
        end = pd.Timestamp(date_str) + pd.Timedelta(days=1)
        end = end.strftime('%Y%m%d%H%M%S')
        ax.axvspan(start, end, alpha=0.2, color='gray')

        # get the daily mean value for the shaded region
        daily_mean_val = daily_mean[date_str]
        daily_mean.index = pd.to_datetime(daily_mean.index)

        arrow_x = pd.to_datetime(start) + pd.to_timedelta((pd.to_datetime(end) - pd.to_datetime(start)) / 2)
        arrow_x = mdates.date2num(arrow_x)  # Convert datetime object to matplotlib date number

        print('date_str:', date_str)
        print('start:', start)
        print('end:', end)
        print('daily_mean_val:', daily_mean_val)
        print('arrow_x:', arrow_x)

        try:
            offset_box = offsetbox.AnnotationBbox(offsetbox.TextArea(f'mean = {daily_mean_val:.2f}'), 
                                                    xy=(arrow_x, daily_mean_val), 
                                                    xybox=(arrow_x, daily_mean_val), 
                                                    boxcoords="data", 
                                                    arrowprops=dict(facecolor='black', arrowstyle='->', 
                                                                    connectionstyle='arc3,rad=0.2', relpos=(0.5, 0.5), 
                                                                    mutation_scale=15), 
                                                    pad=0.5)
            ax.add_artist(offset_box)
        except Exception as e:
            print('Error adding annotation:', e)
file_path = '/ocean/projects/atm200005p/ding0928/nc_file_full/Storm_Peak_julaug_2014.nas'
df = read_data(file_path)
ax = plot_data(df)
add_shaded_regions(ax, df)
plt.show()
plt.savefig('output_fig/obser_storm_peak.png', dpi=800)


