clc;close all;clear;
%% Setting up kernel matrix
rng(3)
n = 150;
m = 200;
TE = linspace(0.3,400,n);
T2 = linspace(1,200,m)';
A=zeros(n,m);
dT = T2(2) - T2(1);

for i=1:n
    for j=1:m
        A(i,j)=exp(-TE(i)/T2(j))*dT; % set up Kernel matrix
    end
end
%% define an arbitrary PDF for testing
npeaks = 2;
T2mu_sim = [30;120];
% T2mu_sim = [30;50];
sigma_i = [3;5];
p = zeros(npeaks,m);
for ii = 1:npeaks
    p(ii,:) = normpdf(T2,T2mu_sim(ii),sigma_i(ii));
end
f_coef = [1/3; 2/3];
IdealModel_weighted = p'*(f_coef);
dat_noiseless = A*IdealModel_weighted;
SNR = 500;
Lambda = logspace(-6,1,16);
load lambda_16_SNR_500_nrun_20_sigma_min_2_sigma_max_4_basis2_1604020lmbda_min-6lmbda_max1.mat
%% run repetitive simulations
nrun = 10;
F_rec_vec = zeros(m,nrun);
DP_rec_vec = zeros(m,nrun);

for i = 1:nrun
    dat_noisy = dat_noiseless + max(abs(dat_noiseless))/SNR*randn(length(TE),1);
    IdealModel_weighted = IdealModel_weighted/max(abs(dat_noisy));
    dat_noisy = dat_noisy/max(abs(dat_noisy));
    %% DP
    [DP_rec_vec(:,i),id_dp] = discrep_L2(dat_noisy,A,SNR,Lambda);
    id_dp
    %% online computation
    [F_rec_vec(:,i),alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum1(dat_noisy,Gaus_info);
end

%% generate figures
figure;
for i = 1:nrun
subplot(3,nrun,i)
plot(T2,IdealModel_weighted,'LineWidth',2,'color',[0,0,0])
hold on
plot(T2,DP_rec_vec(:,i),'LineWidth',2,'color',[0, 0, 1]);
hold off
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)

subplot(3,nrun,nrun+i)
plot(T2,IdealModel_weighted,'LineWidth',2,'color',[0,0,0])
hold on
plot(T2,F_rec_vec(:,i),'LineWidth',2,'color',[1, 0, 0]);
hold off
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)

subplot(3,nrun,2*nrun+i)
plot(TE,A*IdealModel_weighted,'LineWidth',2,'color',[0,0,0])
hold on
plot(TE,A*DP_rec_vec(:,i),':','LineWidth',2,'color',[0, 0, 1]);
plot(TE,A*F_rec_vec(:,i),'--','LineWidth',2,'color',[1, 0, 0]);
hold off
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
end
set(gcf,'position',[32         154        1824         779]);
% saveas(gcf,'figures/resolve_peaks_comparison','epsc')
% saveas(gcf,'figures/resolve_peaks_comparison','png')
saveas(gcf,'figures/resolve_peaks_comparison_well_posed','epsc')
saveas(gcf,'figures/resolve_peaks_comparison_well_posed','png')