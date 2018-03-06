% by Shu Wu: swu33@wisc.edu
% 2/1/2016
addpath('/Users/shu/work/peru/deliverable/codes/PR')
addpath('/Users/shu/work/peru/deliverable/codes/HGT')
addpath('/Users/shu/work/peru/deliverable/codes/TP')
addpath('/Users/shu/work/peru/deliverable/data')
clear;clc;close all;
N_yr = 49;
% Precip
load PR_EOF_12.mat
T1(1:2,:) = T(1:2,1:N_yr-1);
load PR_EOF_1.mat
T2(1:2,:) = T(1:2,2:N_yr);
load PR_EOF_2.mat
T3(1:2,:) = T(1:2,2:N_yr);
load PR_EOF_3.mat
T4(1:2,:) = T(1:2,2:N_yr);
% HGT 200 
load NCEP1_HGT_GLB200_12.mat
T1(3:4,:) = T(1:2,1:N_yr-1);
load NCEP1_HGT_GLB200_1.mat
T2(3:4,:) = T(1:2,2:N_yr);
load NCEP1_HGT_GLB200_2.mat
T3(3:4,:) = T(1:2,2:N_yr);
load NCEP1_HGT_GLB200_3.mat
T4(3:4,:) = T(1:2,2:N_yr);
% SST starts from 1965
load TP_EOF_mn_12.mat	
T1(5:7,:) = T(1:3,2:N_yr);
load TP_EOF_mn_1.mat
T2(5:7,:) = T(1:3,3:N_yr+1);
load TP_EOF_mn_2.mat
T3(5:7,:) = T(1:3,3:N_yr+1);
load TP_EOF_mn_3.mat
T4(5:7,:) = T(1:3,3:N_yr+1);
N_yr = N_yr - 1;

X = cat(1,T1,T2,T3,T4);

SNS = 4;
[B,G,N] = SLIM(X,SNS);

% forecast
T1_P = squeeze(G(:,:,1))*T1;
T1_P = T2;
T2_P = squeeze(G(:,:,2))*T1_P;
T3_P = squeeze(G(:,:,3))*T2_P;

% load precipitation
load PCP_ALL_1966_2016.mat
PCP = squeeze(nanmean(PCP_ALL,3));  
SPCCidx = squeeze(nansum(PCP(1:3,:)));
SPCCidx_stn = squeeze(nansum(PCP_ALL(1:3,:,:)));

mn = 2;
%SPCCidx = squeeze((PCP(mn,:)));
%SPCCidx_stn = squeeze((PCP_ALL(mn,:,:)));
%TR = squeeze(nanmean(PCP(1:3,2:end)));  % ! right
%PR = squeeze(PCP(2,2:end));
yr_b = 1966;
yr_e = 2015;

load PR_EOF_1.mat
PR_JAN_HD = L(:,1:2)*T1_P(1:2,:);

load PR_EOF_2.mat
PR_FEB_HD = L(:,1:2)*T2_P(1:2,:);

load PR_EOF_3.mat
PR_MAR_HD = L(:,1:2)*T3_P(1:2,:);

MEANP = (PR_JAN_HD+PR_FEB_HD+PR_MAR_HD);


%MEANP = PR_FEB_HD;

SPCCidx_p = squeeze(mean(MEANP));
MEAN = mean(SPCCidx);
STD = std(SPCCidx);

for im = 1:29
    temp = corrcoef(MEANP(im,1:end-1),SPCCidx_stn(2:length(SPCCidx_p),im)');
    R(im) = temp(1,2);
end


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
yr = 1967:2016;
subplot(2,1,1)
plot(yr,SPCCidx(2:end),'k','linewidth',3)
hold on
plot(yr(1:length(SPCCidx_p)),SPCCidx_p+MEAN,'--k','linewidth',3)
legend('Observed SPPidx','Predicted SPPidx');
plot(yr,ones(1,length(yr))*(MEAN),':k','linewidth',3)
plot(yr,ones(1,length(yr))*(MEAN+STD),':k','linewidth',3)
plot(yr,ones(1,length(yr))*(MEAN-STD),':k','linewidth',3)
set(gca,'fontsize',30)
set(gcf,'color','w')
xlim([yr(1) yr(end)])
temp = corrcoef(SPCCidx(2:length(SPCCidx_p)),SPCCidx_p(1:end-1));
tmp = (SPCCidx(2:length(SPCCidx_p))-SPCCidx_p(1:end-1));
RMSE = sqrt(mean(tmp.^2))
ylabel('SPPidx(mm)')
title(sprintf('COR=%2.2f',temp(1,2)))
set(gcf,'paperpositionmode','auto')
%print('LIM_SKILL','-dpng');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*9/10])
bar(R)
set(gca,'xtick',1:29)
set(gca,'fontsize',30)
set(gcf,'color','w')
