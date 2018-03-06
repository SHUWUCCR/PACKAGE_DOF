% test_tmp

clear;clc;

file = '/Users/shu/data/HadISST_sst.nc.gz.1'

file = '/Users/shu/data/CenTrends/CenTrends_v1_monthly.nc';


Ncinfo_DOF(file);


target_varName = 'precip'
target_timeName = 'time'
target_lonName = 'longitude'
target_latName = 'latitude'
target_lon = [34 40];
target_lat = [6 14];


[output lon lat] = Ncread_DOF(file, ...
target_varName,target_timeName,target_lonName,target_latName,target_lon,target_lat);


