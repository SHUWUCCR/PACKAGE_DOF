clear;clc;tic;

%1870-now
file_sst =  '/Users/shu/data/HADISST/HadISST_sst.nc';
target_varName = 'sst';
target_timeName = 'time';
target_lonName = 'longitude';
target_latName = 'latitude';

%Tropical Pacific
if 1
   target_lon = [120 270];
   target_lat = [-20 20];
end

[sst lon_sst lat_sst] = Ncread_DOF(file_sst, ...
    target_varName,target_timeName,target_lonName,target_latName,target_lon,target_lat);

yr_b = 1981;
yr_e = 2016;

sst = sst(:,:,(yr_b-1870)*12+1:end);
sst_AMJJAS = zeros(size(sst,1),size(sst,2),6,nyr);
for iyr = 1:nyr
    sst_AMJJAS(:,:,:,iyr) = sst(:,:,(iyr-1)*12+4:(iyr-1)*12+9);
end
clear sst
%%%% end of reading data%%%%%%%%%%%

% caculate anomaly
mean_sst = squeeze(mean(sst_AMJJAS,4));
for iyr = 1:nyr
    sst_AMJJAS(:,:,:,iyr) = squeeze(sst_AMJJAS(:,:,:,iyr)) - mean_sst;
end
% add weight
PI = 3.1415926;
wgt_sst = cos(lat_sst*PI/180);
sst_AMJJAS_bk = sst_AMJJAS;
for im = 1:length(lat_sst)
    sst_AMJJAS(:,im,:,:) = sst_AMJJAS(:,im,:,:)*wgt_sst(im);
end

% calculate EOF
n_lon_sst = length(lon_sst);
n_lat_sst = length(lat_sst);
sst_AMJJAS = reshape(sst_AMJJAS,n_lon_sst*n_lat_sst,6*nyr);
sst_AMJJAS_bk = reshape(sst_AMJJAS_bk,n_lon_sst*n_lat_sst,6*nyr);
missing_sst = isnan(squeeze(mean(sst_AMJJAS,2)));
sst_AMJJAS(missing_sst,:) = [];
sst_AMJJAS_bk(missing_sst,:) = [];
[contri_sst,l_sst] = EOF(sst_AMJJAS);
t_sst = l_sst'*sst_AMJJAS_bk;

toc;



