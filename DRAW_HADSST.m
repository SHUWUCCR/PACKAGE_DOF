% by Shu Wu: swu33@wisc.edu
% 11/15/2016

clear;clc;close all

ncid = netcdf.open('/Users/shu/data/HadISST_sst.nc.gz','nowrite');
varid = netcdf.inqVarID(ncid,'sst');
data = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'longitude');
lon = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'latitude');
lat = netcdf.getVar(ncid,varid);
netcdf.close(ncid);
data = double(data);
lat = double(lat);
lon = double(lon);

% add path of figure tools
addpath('/Users/shu/m_map')

HDyr_b = 1870;
yr_b = 1900;
yr_e = 2014;
n_yr = yr_e - yr_b + 1;
yri = 1982;
mni = 5;

ssti = squeeze(data(:,:,(yri-HDyr_b)*12+mni));
data = data(:,:,(yr_b-HDyr_b)*12+1:(yr_e-HDyr_b+1)*12);
MEAN = squeeze(mean(data(:,:,mni:12:end),3));

ssti(abs(ssti)>1e3) = NaN;
MEAN(abs(MEAN)>1e3) = NaN';

ssta = ssti - MEAN;

lat_b = 60;
lat_e = -60;
[C lat_B] = min(abs(lat_b-lat));
[C lat_E] = min(abs(lat_e-lat));
lat = lat(lat_B:lat_E);
ssta = ssta(:,lat_B:lat_E);
[C lon_B] = min(abs(lon-30));
lon = lon(lon_B):lon(lon_B)+359;
ssta2 = cat(1,ssta(lon_B:end,:),ssta(1:lon_B-1,:));
ssta = fliplr(ssta2);
lat = flipud(lat(:));
 
    ll =  linspecer;
    ca = 4;
    inv= 0.5;
    DL = -ca:inv:ca;
    INV = floor(size(ll,1)/length(DL));
    lineStyles = (ll(1:INV:end,:));
    lineStyles =  lineStyles(1:length(DL)-1,:);
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
        temp = ssta;
        m_proj('miller','long',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
        hold on
        colormap(lineStyles)
        [cc,hh]=m_contourf(lon,lat,temp',DL);
        set(hh,'lineStyle','none');
        [cc,hh]=m_contour(lon,lat,temp',[0 0],'k','linewidth',3);
        m_coast('line','color',[0.1 0.1 0.1]);
        hcb=colorbar;
        set(hcb,'YTick',[DL])
        set(gcf,'color','w')
        set(gcf,'paperpositionmode','auto')
        caxis([-ca ca])
        set(gca,'fontsize',30)
        m_grid('linestyle','none','tickdir','out','linewidth',3);
  
    set(gcf,'paperpositionmode','auto')
    


