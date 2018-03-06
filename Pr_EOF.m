clear;clc;tic;close all;

%1981-2017
file_pr =  '/Users/shu/data/CHIRPS/new/chirps-v2.0.monthly.nc';
target_varName = 'precip';
target_timeName = 'time';
target_lonName = 'longitude';
target_latName = 'latitude';

%Shukla Ethiopia
if 1
    NAME = 'PR_SETH'
    NAME2 = [NAME '_AMJJAS'];
    NAME3 = [NAME '_AMJJAS.jpg'];
    target_lat = [6 14];
    target_lon = [34 40];
end
[pr pr_lon pr_lat] = Ncread_DOF(file_pr, ...
    target_varName,target_timeName,target_lonName,target_latName,target_lon,target_lat);
% cut the data for AMJJAS from 1981-2017
nyr = floor(size(pr,3)/12);
% get AMJJAS data
pr_AMJJAS = zeros(size(pr,1),size(pr,2),6,nyr);

for iyr = 1:nyr
    pr_AMJJAS(:,:,:,iyr) = pr(:,:,(iyr-1)*12+4:(iyr-1)*12+9);
end
clear pr

% caculate anomaly
mean_pr = squeeze(mean(pr_AMJJAS,4));
for iyr = 1:nyr
    pr_AMJJAS(:,:,:,iyr) = squeeze(pr_AMJJAS(:,:,:,iyr)) - mean_pr;
end
% add weight
PI = 3.1415926;
wgt_pr = cos(pr_lat*PI/180);
pr_AMJJAS_bk = pr_AMJJAS;
for im = 1:length(pr_lat)
    pr_AMJJAS(:,im,:,:) = pr_AMJJAS(:,im,:,:)*wgt_pr(im);
end
% calculate EOF
n_pr_lon = length(pr_lon);
n_pr_lat = length(pr_lat);
pr_AMJJAS = reshape(pr_AMJJAS,n_pr_lon*n_pr_lat,6*nyr);
pr_AMJJAS_bk = reshape(pr_AMJJAS_bk,n_pr_lon*n_pr_lat,6*nyr);
missing_pr = isnan(squeeze(mean(pr_AMJJAS,2)));
pr_AMJJAS(missing_pr,:) = [];
pr_AMJJAS_bk(missing_pr,:) = [];
[pr_contri,pr_l] = EOF(pr_AMJJAS);
pr_t = pr_l'*pr_AMJJAS_bk;

addpath('/Users/shu/m_map')
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
ll =  linspecer;
ca = 60;
inv = 10;
DL = -ca:inv:ca;
INV = floor(size(ll,1)/length(DL));
lineStyles = (ll(1:INV:end,:));
lineStyles =  lineStyles(1:length(DL)-1,:);

save(NAME2,'pr_contri','pr_l','pr_t', ...
    'pr_AMJJAS_bk','pr_lon','pr_lat');

for im = 1:3
    subplot(1,3,im)
    temp = pr_l(:,im)*std(pr_t(im,:),1);
    temp = reshape(temp,n_pr_lon,n_pr_lat)';
    m_proj('miller','long',[pr_lon(1) pr_lon(end)],'lat',[pr_lat(1) pr_lat(end)]);
    hold on
    colormap(lineStyles)
    [cc,hh]=m_contourf(pr_lon,pr_lat,temp,DL);
    set(hh,'lineStyle','none');
    [cc,hh]=m_contour(pr_lon,pr_lat,temp,[0 0],'k','linewidth',3);
    m_coast('line','color',[0.1 0.1 0.1]);
    titlename=strcat(num2str(round(pr_contri(im)*10000)/100),'% variance explained');
    title(titlename,'fontsize',20);
    hcb=colorbar;
    set(hcb,'YTick',[DL])
    set(gcf,'color','w')
    set(gcf,'paperpositionmode','auto')
    caxis([-ca ca])
    set(gca,'fontsize',12)
    m_grid('linestyle','none','tickdir','out','linewidth',3);
end
set(gcf,'paperpositionmode','auto')
saveas(gcf,NAME3,'jpg')
%print(NAME3,'jpeg')


toc;



