% by Shu Wu: swu33@wisc.edu
% 2/1/2016

clear;clc;close all;
addpath('/Users/shu/m_map')
N_yr = 36;
% Precip start from 1981

mi = 1:2;
ti = 1:2;

% start from 1981
load SETH_5
T1(mi,:) = T(ti,1:N_yr);
load SETH_6
T2(mi,:) = T(ti,1:N_yr);
load SETH_7
T3(mi,:) = T(ti,1:N_yr);
load SETH_8
T4(mi,:) = T(ti,1:N_yr);
load SETH_9
T5(mi,:) = T(ti,1:N_yr);



% SST starts from 1965
if 1
    % tropical Pacific
    t_b = 1981-1965+1;
    mi = 3;
    ti = 1;
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_5.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_6.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_7.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_8.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_9.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
end

% add NP
% SST starts from 1965
t_b = 1981-1965+1;
mi = 4;
ti = 1;
load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_5.mat
T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_6.mat
T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_7.mat
T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_8.mat
T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_9.mat
T5(mi,:) = T(ti,t_b:t_b+N_yr-1);

if 1
    % add TA
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 5;
    ti = 1;
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_5.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_6.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_7.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_8.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_9.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
end

if 1
    % add NA
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 6;
    ti = 1;
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_5.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_6.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_7.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_8.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_9.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
end


if 1
    % add IN
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 7;
    ti = 1;
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_5.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_6.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_7.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_8.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_9.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
end


if 1
    % HGT 200 start from 1966
    t_b = 1981-1966+1;
    mi = 8;
    ti = 1;
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_5.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_6.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_7.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_8.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_9.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    
end

if 1
    % HGT 850 start from 1966
    t_b = 1981-1966+1;
    mi = 9;
    ti = 1;
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_5.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_6.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_7.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_8.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_9.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
end

if 0
    N_yr = size(T1,2);
    T_tau = T2*T1'/N_yr;
    T0 = T1*T1'/N_yr;
    G_tau1 = T_tau*inv(T0);
    B_tau1 = real(logm(G_tau1));
    N1 = -B_tau1*T0+T0*B_tau1;
    T_tau = T3*T2'/N_yr;
    T0 = T2*T2'/N_yr;
    G_tau2 = T_tau*inv(T0);
    B_tau2 = real(logm(G_tau2));
    N2 = -B_tau2*T0+T0*B_tau2;
    T_tau = T4*T3'/N_yr;
    T0 = T3*T3'/N_yr;
    G_tau3 = T_tau*inv(T0);
    B_tau3 = real(logm(G_tau3));
    N3 = -B_tau3*T0+T0*B_tau3;
    T_tau = T5*T4'/N_yr;
    T0 = T4*T4'/N_yr;
    G_tau4 = T_tau*inv(T0);
    B_tau4 = real(logm(G_tau3));
    N4 = -B_tau4*T0+T0*B_tau4;
end

X = cat(1,T1,T2,T3,T4,T5);

SNS = 5;

[B,G,N] = SLIM(X,SNS);

% forecast
T1_P = squeeze(G(:,:,1))*T1;
T2_P = squeeze(G(:,:,2))*T1_P;
T3_P = squeeze(G(:,:,3))*T2_P;
T4_P = squeeze(G(:,:,4))*T3_P;


if 0
    T1_P = expm(squeeze(B(:,:,1)))*T1;
    T2_P = expm(squeeze(B(:,:,2)))*T1_P;
    T3_P = expm(squeeze(B(:,:,3)))*T2_P;
    T4_P = expm(squeeze(B(:,:,4)))*T3_P;
    
end

if 0
    T1_P = G_tau1*T1;
    T2_P = G_tau2*T1_P;
    T3_P = G_tau3*T2_P;
    T4_P = G_tau4*T3_P;
    
    % real part of VB needs to be negative
    [LB1 VB1] = eig(B_tau1);
    [LB2 VB2] = eig(B_tau2);
    [LB3 VB3] = eig(B_tau3);
    [LB4 VB4] = eig(B_tau4);
    
    [LN1 VN1] = eig(N1);
    [LN2 VN2] = eig(N2);
    [LN3 VN3] = eig(N3);
    [LN4 VN4] = eig(N4);
    
    
end
%G_tau1 = expm(B_tau1);
%G_tau2 = expm(B_tau2);
%G_tau3 = expm(B_tau3);

%corrcoef(T1_P(1,:),T2(1,:))
%corrcoef(T2_P(1,:),T3(1,:))
%corrcoef(T3_P(1,:),T4(1,:))


% load precipitation
load PR_SETH
PR = reshape(PR,N_LON*N_LAT,12,N_yr);
PCP = squeeze(nanmean(PR,1));
SPCCidx = squeeze(nansum(PCP(6:9,:)));
%TR = squeeze(nanmean(PCP(1:3,2:end)));  % ! right

yr_b = 1981;
yr_e = 2016;

mi = 1:2;
load SETH_6
PR_JUN_HD = L(:,mi)*T1_P(mi,:);

load SETH_7
PR_JUL_HD = L(:,mi)*T2_P(mi,:);

load SETH_8
PR_AUG_HD = L(:,mi)*T3_P(mi,:);

load SETH_9
PR_SEP_HD = L(:,mi)*T4_P(mi,:);


MEANP = (PR_JUN_HD+PR_JUL_HD+PR_AUG_HD+PR_SEP_HD);
SPCCidx_p = squeeze(mean(MEANP));
MEAN = mean(SPCCidx);
STD = std(SPCCidx);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
yr = 1981:2016;

plot(yr,SPCCidx(1:end),'k','linewidth',3)
hold on
plot(yr(1:length(SPCCidx_p)),SPCCidx_p+MEAN,'g','linewidth',3)
legend('Observed ETHidx','Predicted ETHidx');
scatter(yr,SPCCidx(1:end),200,'filled')
plot(yr,ones(1,length(yr))*(MEAN),'--r','linewidth',2)
plot(yr,ones(1,length(yr))*(MEAN+STD),'--r','linewidth',2)
plot(yr,ones(1,length(yr))*(MEAN-STD),'--r','linewidth',2)
set(gca,'fontsize',30)
set(gcf,'color','w')
set(gca,'xtick',1981:5:2017)
temp = corrcoef(SPCCidx,SPCCidx_p);
ylabel('ETHidx')
title(sprintf('COR=%2.2f',temp(1,2)))

%
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
ll =  linspecer;
ca = 1;
inv = 0.2;
DL = -ca:inv:ca;
INV = floor(size(ll,1)/length(DL));
lineStyles = (ll(1:INV:end,:));
lineStyles =  lineStyles(1:length(DL)-1,:);


for m = 1:size(PR,1)
    temp = corrcoef(squeeze(PR(m,6,:)),PR_JUN_HD(m,:));
    R6(m) = temp(1,2);
    temp = corrcoef(squeeze(PR(m,7,:)),PR_JUL_HD(m,:));
    R7(m) = temp(1,2);
    temp = corrcoef(squeeze(PR(m,8,:)),PR_AUG_HD(m,:));
    R8(m) = temp(1,2);
    temp = corrcoef(squeeze(PR(m,9,:)),PR_SEP_HD(m,:));
    R9(m) = temp(1,2);
end



subplot(2,2,1)
temp = R6;
temp = reshape(temp,N_LON,N_LAT)';
m_proj('miller','long',[lon(1) lon(length(lon(:)))],'lat',[lat(1) lat(length(lat(:)))]);
hold on
colormap(lineStyles)
[cc,hh]=m_contourf(lon,lat,temp,DL);
set(hh,'lineStyle','none');
[cc,hh]=m_contour(lon,lat,temp,[0 0],'k','linewidth',3);
m_coast('line','color',[0.1 0.1 0.1]);
hcb=colorbar;
set(hcb,'YTick',[DL])
set(gcf,'color','w')
set(gcf,'paperpositionmode','auto')
caxis([-ca ca])
set(gca,'fontsize',12)
m_grid('linestyle','none','tickdir','out','linewidth',3);

subplot(2,2,2)
temp = R7;
temp = reshape(temp,N_LON,N_LAT)';
m_proj('miller','long',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
hold on
colormap(lineStyles)
[cc,hh]=m_contourf(lon,lat,temp,DL);
set(hh,'lineStyle','none');
[cc,hh]=m_contour(lon,lat,temp,[0 0],'k','linewidth',3);
m_coast('line','color',[0.1 0.1 0.1]);
hcb=colorbar;
set(hcb,'YTick',[DL])
set(gcf,'color','w')
set(gcf,'paperpositionmode','auto')
caxis([-ca ca])
set(gca,'fontsize',12)
m_grid('linestyle','none','tickdir','out','linewidth',3);

subplot(2,2,3)
temp = R7;
temp = reshape(temp,N_LON,N_LAT)';
m_proj('miller','long',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
hold on
colormap(lineStyles)
[cc,hh]=m_contourf(lon,lat,temp,DL);
set(hh,'lineStyle','none');
[cc,hh]=m_contour(lon,lat,temp,[0 0],'k','linewidth',3);
m_coast('line','color',[0.1 0.1 0.1]);
hcb=colorbar;
set(hcb,'YTick',[DL])
set(gcf,'color','w')
set(gcf,'paperpositionmode','auto')
caxis([-ca ca])
set(gca,'fontsize',12)
m_grid('linestyle','none','tickdir','out','linewidth',3);

subplot(2,2,4)
temp = R9;
temp = reshape(temp,N_LON,N_LAT)';
m_proj('miller','long',[lon(1) lon(end)],'lat',[lat(1) lat(end)]);
hold on
colormap(lineStyles)
[cc,hh]=m_contourf(lon,lat,temp,DL);
set(hh,'lineStyle','none');
[cc,hh]=m_contour(lon,lat,temp,[0 0],'k','linewidth',3);
m_coast('line','color',[0.1 0.1 0.1]);
hcb=colorbar;
set(hcb,'YTick',[DL])
set(gcf,'color','w')
set(gcf,'paperpositionmode','auto')
caxis([-ca ca])
set(gca,'fontsize',12)
m_grid('linestyle','none','tickdir','out','linewidth',3);





%      print(NAME_save3,'-dpng')
set(gcf,'paperpositionmode','auto')





