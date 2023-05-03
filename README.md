the main branch is after I moved everthing to the ocean directory, and tested successfully with the subregion domain;
Now I want to modify the figure to be 1*2(left: high altitude; right: low altitude)

# 3 May update
obs_groudn_storm_peak.py is to analyze the observed data for the whole time series;
time_series_dask_working.py  time_series_dask_working.ipynb is the long python script for the model output, which is still very slow;

# working logic: still under optimization 
time_series_dask_working.ipynb is where I write the code for convenence reason, 
but interact -n 32c is much faster to run time_series_dask_working.py in terminal, 
so pair terminal with vscode jupyter is by far the fast way to test the code;
