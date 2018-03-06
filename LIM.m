clear;clc;tic;

%1981-2017
file_pr =  '/Users/shu/data/CHIRPS/new/chirps-v2.0.monthly.nc';
target_varName = 'precip'
target_timeName = 'time'
target_lonName = 'longitude'
target_latName = 'latitude'
%Shukla Ethiopia
if 1
   target_lat = [6 14];
   target_lon = [34 40];
end
[pr lon_pr lat_pr] = Ncread_DOF(file_pr, ...
    target_varName,target_timeName,target_lonName,target_latName,target_lon,target_lat);
% cut the data for AMJJAS from 1981-2017
nyr = floor(size(pr,3)/12);
% get AMJJAS data
pr_AMJJAS = zeros(size(pr,1),size(pr,2),6,nyr);

for iyr = 1:nyr
    pr_AMJJAS(:,:,:,iyr) = pr(:,:,(iyr-1)*12+4:(iyr-1)*12+9);
end
clear pr
%1870-now
file_sst =  '/Users/shu/data/HADISST/HadISST_sst.nc';
target_varName = 'sst'
target_timeName = 'time'
target_lonName = 'longitude'
target_latName = 'latitude'
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
mean_pr = squeeze(mean(pr_AMJJAS,4));
mean_sst = squeeze(mean(sst_AMJJAS,4));
for iyr = 1:nyr
    pr_AMJJAS(:,:,:,iyr) = squeeze(pr_AMJJAS(:,:,:,iyr)) - mean_pr;
    sst_AMJJAS(:,:,:,iyr) = squeeze(sst_AMJJAS(:,:,:,iyr)) - mean_sst;
end
% add weight
PI = 3.1415926;
wgt_pr = cos(lat_pr*PI/180);
wgt_sst = cos(lat_sst*PI/180);
pr_AMJJAS_bk = pr_AMJJAS;
sst_AMJJAS_bk = sst_AMJJAS;
for im = 1:length(lat_pr)
    pr_AMJJAS(:,im,:,:) = pr_AMJJAS(:,im,:,:)*wgt_pr(im);
end
for im = 1:length(lat_sst)
    sst_AMJJAS(:,im,:,:) = sst_AMJJAS(:,im,:,:)*wgt_sst(im);
end
% calculate EOF
n_lon_pr = length(lon_pr);
n_lat_pr = length(lat_pr);
pr_AMJJAS = reshape(pr_AMJJAS,n_lon_pr*n_lat_pr,6*nyr);
pr_AMJJAS_bk = reshape(pr_AMJJAS_bk,n_lon_pr*n_lat_pr,6*nyr);
missing_pr = isnan(squeeze(mean(pr_AMJJAS,2)));
pr_AMJJAS(missing_pr,:) = [];
pr_AMJJAS_bk(missing_pr,:) = [];
[contri_pr,l_pr] = EOF(pr_AMJJAS);
t_pr = l_pr'*pr_AMJJAS_bk;

n_lon_sst = length(lon_sst);
n_lat_sst = length(lat_sst);
sst_AMJJAS = reshape(sst_AMJJAS,n_lon_sst*n_lat_sst,6*nyr);
sst_AMJJAS_bk = reshape(sst_AMJJAS_bk,n_lon_sst*n_lat_sst,6*nyr);
missing_sst = isnan(squeeze(mean(sst_AMJJAS,2)));
sst_AMJJAS(missing_sst,:) = [];
sst_AMJJAS_bk(missing_sst,:) = [];
[contri_sst,l_sst] = EOF(sst_AMJJAS);
t_sst = l_sst'*sst_AMJJAS_bk;
% for LIM
x = [t_pr(1:3,:);t_sst(1:3,:)];
x = reshape(x,6,6,nyr);
x0 = x(:,1:5,:);
x1 = x(:,2:end,:);
x0 = reshape(x0,6,5*nyr);
x1 = reshape(x1,6,5*nyr);
g = x1*x0'*inv(x0*x0');
x1_p = g*x0;

for im = 1:6
    temp = corrcoef(x1_p(im,:),x1(im,:));
    cor(im) = temp(1,2);
end

bar(cor)

toc;



