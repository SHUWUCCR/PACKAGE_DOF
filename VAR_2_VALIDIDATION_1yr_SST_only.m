% by Shu Wu: swu33@wisc.edu
% 2/1/2016

clear;clc;close all;
addpath('/Users/shu/m_map')
% N_yr = 20;
% N_yr = 36;
N_yr = 36;
% Precip start from 1981

load PR_SETH
PR = reshape(PR,N_LON*N_LAT,12,N_yr);
PCP = squeeze(nanmean(PR,1));
PCP_BK = PCP;


if 1

mi = 1:2;
ti = 1:2;

% start from 1981

load SETH_6
tpr(mi,:,1) = T(ti,1:N_yr);
load SETH_7
tpr(mi,:,2) = T(ti,1:N_yr);
load SETH_8
tpr(mi,:,3) = T(ti,1:N_yr);
load SETH_9
tpr(mi,:,4) = T(ti,1:N_yr);
tpr_bk = tpr;
end


if 1
    % SST starts from 1965
    % tropical Pacific
    t_b = 1981-1965+1;
    mi = 1:2;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TP/TP_EOF_mn_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end


if 1
    % add TA
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 3:4;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/TA/TA_EOF_mn_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end


if 1
    % add IN
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 5:6;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/IN/IN_EOF_mn_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end


if 0
    % HGT 200 start from 1966
    t_b = 1981-1966+1;
    mi = 7:8;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB200_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end


if 0
    % HGT 850 start from 1966
    t_b = 1981-1966+1;
    mi = 11:12;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/HGT/NCEP1_HGT_GLB850_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end

if 1
    % add NP
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 7:8;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NP/NP_EOF_mn_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end

if 1
    % add NA
    % SST starts from 1965
    t_b = 1981-1965+1;
    mi = 9:10;
    mi_e = 10;
    ti = 1:2;
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_4.mat
    T1(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_5.mat
    T2(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_6.mat
    T3(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_7.mat
    T4(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_8.mat
    T5(mi,:) = T(ti,t_b:t_b+N_yr-1);
    load /Users/shu/work/peru/deliverable/codes/NA/NA_EOF_mn_9.mat
    T6(mi,:) = T(ti,t_b:t_b+N_yr-1);
end

if 0
T1([12],:) = [];
T2([12],:) = [];
T3([12],:) = [];
T4([12],:) = [];
T5([12],:) = [];
T6([12],:) = [];
end

T12 = cat(1,T1,T2);
T23 = cat(1,T2,T3);
T34 = cat(1,T3,T4);
T45 = cat(1,T4,T5);
T56 = cat(1,T5,T6);

T1_BK = T12;
T2_BK = T23;
T3_BK = T34;
T4_BK = T45;
T5_BK = T56;

N_yr_BK = N_yr;
TP1 = zeros(mi_e,N_yr);
TP2 = zeros(mi_e,N_yr);
TP3 = zeros(mi_e,N_yr);
TP4 = zeros(mi_e,N_yr);
% cross valide test start
for LF = 1:36;
    LF
    T1 = T1_BK;
    T2 = T2_BK;
    T3 = T3_BK;
    T4 = T4_BK;
    T5 = T5_BK;
    T00 = T1(:,LF);
    tpr = tpr_bk;
    tpr(:,LF,:) = [];
    
    PCP = PCP_BK;
    PCP(:,LF) = [];
    

    T1(:,LF) = [];
    T2(:,LF) = [];
    T3(:,LF) = [];
    T4(:,LF) = [];
    T5(:,LF) = [];
    
    
    for im = 1:2
    xmr = T3(1:mi_e,:);
    ymr = squeeze(tpr(im,:,1));
    b(:,im,1) = regress(ymr',xmr');    
    xmr = T4(1:mi_e,:);
    ymr = squeeze(tpr(im,:,2));
    b(:,im,2) = regress(ymr',xmr');
     xmr = T5(1:mi_e,:);
    ymr = squeeze(tpr(im,:,3));
    b(:,im,3) = regress(ymr',xmr');  
     xmr = T5(mi_e+1:end,:);
    ymr = squeeze(tpr(im,:,4));
    b(:,im,4) = regress(ymr',xmr'); 
    end
    
    xmr = T3(1:mi_e,:);
    ymr = squeeze(PCP(6,:));
    bb(:,1) = regress(ymr',xmr');
    xmr = T4(1:mi_e,:);
    ymr = squeeze(PCP(7,:));
    bb(:,2) = regress(ymr',xmr');  
    xmr = T5(1:mi_e,:);
    ymr = squeeze(PCP(8,:));
    bb(:,3) = regress(ymr',xmr');
    xmr = T5(mi_e+1:end,:);
    ymr = squeeze(PCP(9,:));
    bb(:,4) = regress(ymr',xmr');
    
    
    if 0
         T00([12 size(T1,1)/2+12]) = [];
T1([12 size(T1,1)/2+12],:) = [];
T2([12 size(T1,1)/2+12],:) = [];
T3([12 size(T1,1)/2+12],:) = [];
T4([12 size(T1,1)/2+12],:) = [];
T5([12 size(T1,1)/2+12],:) = [];
    end

    
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
    
    %save G_tau_VAR2_34yr G_tau1 G_tau2 G_tau3 G_tau4
    %stop
    %load G_tau_VAR2_34yr
    
    T1_P = G_tau1*T00;
    T2_P = G_tau2*T1_P;
    T3_P = G_tau3*T2_P;
    T4_P = G_tau4*T3_P;
    
    
    T1_P = T1_P(size(T1_P,1)/2+1:end,:);
    T2_P = T2_P(size(T2_P,1)/2+1:end,:);
    T3_P = T3_P(size(T3_P,1)/2+1:end,:);
    T4_P = T4_P(size(T4_P,1)/2+1:end,:);
    
    % real part of VB needs to be negative
    %[LB1 VB1] = eig(B_tau1);
    %[LB2 VB2] = eig(B_tau2);
    %[LB3 VB3] = eig(B_tau3);
    %[LB4 VB4] = eig(B_tau4);
    
    %[LN1 VN1] = eig(N1);
    %[LN2 VN2] = eig(N2);
    %[LN3 VN3] = eig(N3);
    %[LN4 VN4] = eig(N4);
    
    
    %G_tau1 = expm(B_tau1);
    %G_tau2 = expm(B_tau2);
    %G_tau3 = expm(B_tau3);
    
    %corrcoef(T1_P(1,:),T2(1,:))
    %corrcoef(T2_P(1,:),T3(1,:))
    %corrcoef(T3_P(1,:),T4(1,:))
    
    if 1
    
    mi = 1:2;
load SETH_6
Pr_P1 = squeeze(b(:,mi,1))'*T1_P;
PR_JUN_HD = L(:,mi)*Pr_P1(mi,:);

load SETH_7
Pr_P2 = squeeze(b(:,mi,2))'*T2_P;
PR_JUL_HD = L(:,mi)*Pr_P2(mi,:);

load SETH_8
Pr_P3 = squeeze(b(:,mi,3))'*T3_P;
PR_AUG_HD = L(:,mi)*Pr_P3(mi,:);

load SETH_9
Pr_P4 = squeeze(b(:,mi,4))'*T4_P;
PR_SEP_HD = L(:,mi)*Pr_P4(mi,:);

MEANP = (PR_JUN_HD+PR_JUL_HD+PR_AUG_HD+PR_SEP_HD);
SPCCidx_p = squeeze(mean(MEANP));

SPCCidx_pp(LF) = SPCCidx_p;

    SPCC_JUN_HD = squeeze(bb(:,1))'*T1_P;
    SPCC_JUL_HD = squeeze(bb(:,2))'*T2_P;
    SPCC_AUG_HD = squeeze(bb(:,3))'*T3_P;
    SPCC_SEP_HD = squeeze(bb(:,4))'*T4_P;
    MEANP = (SPCC_JUN_HD+SPCC_JUL_HD+SPCC_AUG_HD+SPCC_SEP_HD);
   SPCCidx_p = squeeze(mean(MEANP));

SPCCidx_pp2(LF) = SPCCidx_p;

    
    end
    
    TP1(:,LF) = T1_P;
    TP2(:,LF) = T2_P;
    TP3(:,LF) = T3_P;
    TP4(:,LF) = T4_P;
end


for im = 1:mi_e
    temp = corrcoef(TP1(im,:),T3_BK(im,:));
    CORS(im,1) = temp(1,2);
    temp = corrcoef(TP2(im,:),T4_BK(im,:));
    CORS(im,2) = temp(1,2);
    temp = corrcoef(TP3(im,:),T5_BK(im,:));
    CORS(im,3) = temp(1,2);
    temp = corrcoef(TP3(im,:),T5_BK(im+mi_e,:));
    CORS(im,4) = temp(1,2);
end

for im = 1:4
    subplot(2,2,im)
    bar(CORS(:,im));
    set(gca,'fontsize',30)
end
set(gcf,'color','w')


SPCCidx_p = SPCCidx_pp;
%SPCCidx_p = SPCCidx_pp2;
load PR_SETH
PR = reshape(PR,N_LON*N_LAT,12,N_yr);
PCP = squeeze(nanmean(PR,1));
SPCCidx = squeeze(nansum(PCP(6:9,:)));
%TR = squeeze(nanmean(PCP(1:3,2:end)));  % ! right

yr_b = 1981;
yr_e = 2016;


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
temp = corrcoef(SPCCidx(1:length(SPCCidx_p)),SPCCidx_p);
ylabel('ETHidx')
title(sprintf('COR=%2.2f',temp(1,2)))







