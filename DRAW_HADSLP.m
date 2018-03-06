% by Shu Wu: swu33@wisc.edu
% 11/15/2016

clear;clc;close all
load /Users/shu/data/SLP/slp.mat

yri = 2007-yr(1)+1;
mni = 12;

yr_b = 1900;
yr_e = 2014;

mslp = squeeze(slp(:,:,mni,yr_b-yr(1)+1:yr_e-yr(1)+1));

MEAN = squeeze(mean(mslp,3));


temp = (squeeze(slp(:,:,mni,yri))-MEAN)';
            
            

N_lat = length(lat);
N_lon = length(lon);
addpath('/Users/shu/m_map')

ll =  linspecer;
%ca = floor(max(abs(temp(:))));
ca = 10;
inv = ca/5;
DL = -ca:inv:ca;
INV = floor(size(ll,1)/length(DL));
lineStyles = (ll(1:INV:end,:));
lineStyles =  lineStyles(1:length(DL)-1,:);
scrsz = get(0,'ScreenSize');
        
        figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
       % subplot(4,1,1:2)
            %L(:,2) = L(:,2)*(-1);

            m_proj('miller','long',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
            hold on
            colormap(lineStyles)
            [cc,hh]=m_contourf(lon,lat,temp,DL);
            set(hh,'lineStyle','none');
            [cc,hh]=m_contour(lon,lat,temp,[0 0],'k','linewidth',3);
          %  m_coast('line','color',[0.1 0.1 0.1],'linewidth',3);
            m_coast('line','color','r','linewidth',3);
            hcb=colorbar;
            set(hcb,'YTick',[DL])
            set(gcf,'color','w')
            set(gcf,'paperpositionmode','auto')
            caxis([-ca ca])
            set(gca,'fontsize',30)
            m_grid('linestyle','none','tickdir','out','linewidth',3);






