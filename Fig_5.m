clc;close all;clear;
rng(1);
%% Setting up kernel matrix
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
T2mu_sim = [30;80];
sigma_i = [3;5];
p = zeros(npeaks,m);
for ii = 1:npeaks
    p(ii,:) = normpdf(T2,T2mu_sim(ii),sigma_i(ii));
end
f_coef = [2/3; 1/3];
IdealModel_weighted = p'*(f_coef);
dat_noiseless = A*IdealModel_weighted;
SNR = 500;

load lambda_16_SNR_500_nrun_20_sigma_min_2_sigma_max_4_basis2_1604020lmbda_min-6lmbda_max1.mat
Lambda_opt = 6*1e-3;
Lambda = Lambda_opt*2.^(-4:1:5);
%% run repetitive simulations
nrun = 1;
F_rec_vec = zeros(m,nrun);
%%
dat_noisy = dat_noiseless + max(abs(dat_noiseless))/SNR*randn(length(TE),1);
IdealModel_weighted = IdealModel_weighted/max(abs(dat_noisy));
dat_noisy = dat_noisy/max(abs(dat_noisy));

%% Conventional Methods chooses one lambda
[x0,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(1));
[x1,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(2));
[x2,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(3));
[x3,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(4));
[x4,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(5));
[x5,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(6));
[x6,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(7));
[x7,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(8));
[x8,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(9));
[x9,mdl_err,xnorm] = nonnegtik_hnorm(A,dat_noisy,Lambda(10));

res_vec = zeros(length(Lambda),1);
res_vec(1) = norm(IdealModel_weighted - x0);
res_vec(2) = norm(IdealModel_weighted - x1);
res_vec(3) = norm(IdealModel_weighted - x2);
res_vec(4) = norm(IdealModel_weighted - x3);
res_vec(5) = norm(IdealModel_weighted - x4);
res_vec(6) = norm(IdealModel_weighted - x5);
res_vec(7) = norm(IdealModel_weighted - x6);
res_vec(8) = norm(IdealModel_weighted - x7);
res_vec(9) = norm(IdealModel_weighted - x8);
res_vec(10) = norm(IdealModel_weighted - x9);

[~,id_OPT] = min(res_vec);

%% online computation
adj_lambda = 2;
[f_rec,alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum2(dat_noisy,id_OPT, adj_lambda, Gaus_info, Lambda);

%%
figure;
hold on
plot(T2, x2,'LineWidth',1.5)
plot(T2, x3,'LineWidth',1.5)
plot(T2, x4,'LineWidth',1.5)
plot(T2, x5,'LineWidth',1.5)
plot(T2, x6,'LineWidth',1.5)
plot(T2, IdealModel_weighted,'LineWidth',2.5,'Color','k')
legend({
    'f_{\lambda_1}',...
    'f_{\lambda_2}',...
    'f_{\lambda_3}(OPTIMAL)',...
    'f_{\lambda_4}',...
    'f_{\lambda_5}',...
    'True Dist.'},'FontSize',16)
ylim([0 0.2])
set(gcf,'position',[1327         841         673         497])
saveas(gcf,'figures/lambda_shift_stability_A','epsc')
saveas(gcf,'figures/lambda_shift_stability_A','png')


% MultiReg chooses multiple lambdas
Alpha = alpha_L2;
X2 = [x0 x1 x2 x3 x4]*Alpha;
X3 = [x1 x2 x3 x4 x5]*Alpha;
X4 = [x2 x3 x4 x5 x6]*Alpha;
X5 = [x3 x4 x5 x6 x7]*Alpha;
X6 = [x4 x5 x6 x7 x8]*Alpha;
% 
% X2 = [x1 x2 x3]*Alpha;
% X3 = [x2 x3 x4]*Alpha;
% X4 = [x3 x4 x5]*Alpha;
% X5 = [x4 x5 x6]*Alpha;
% X6 = [x5 x6 x7]*Alpha;

figure;
hold on
plot(T2, X2,'LineWidth',1.5)
plot(T2, X3,'LineWidth',1.5)
plot(T2, X4,'LineWidth',1.5)
plot(T2, X5,'LineWidth',1.5)
plot(T2, X6,'LineWidth',1.5)
plot(T2, IdealModel_weighted,'LineWidth',2.5,'Color','k')

legend({
    '\alpha_1 f_{\lambda_0}+ \alpha_2 f_{\lambda_1}+ \alpha_3f_{\lambda_2}',...
    '\alpha_1 f_{\lambda_1}+ \alpha_2 f_{\lambda_2}+ \alpha_3f_{\lambda_3}',...
    '\alpha_1 f_{\lambda_2}+ \alpha_2 f_{\lambda_3}+ \alpha_3f_{\lambda_4} (SpanReg)',...
    '\alpha_1 f_{\lambda_3}+ \alpha_2 f_{\lambda_4}+ \alpha_3f_{\lambda_5}',...
    '\alpha_1 f_{\lambda_4}+ \alpha_2 f_{\lambda_5}+ \alpha_3f_{\lambda_6}',...
    'True Dist.'},'FontSize',16)
ylim([0 0.2])
set(gcf,'position',[1327         841         673         497])
saveas(gcf,'figures/lambda_shift_stability_B','epsc')
saveas(gcf,'figures/lambda_shift_stability_B','png')


%% find differences 
dx1 = x4 - x2;
dx2 = x4 - x3;
dx3 = x4 - x5;
dx4 = x4 - x6;

figure;
plot(T2, dx1,'LineWidth',1.5)
hold on
plot(T2, dx2,'LineWidth',1.5)
plot(T2, dx3,'LineWidth',1.5)
plot(T2, dx4,'LineWidth',1.5)
legend({'\Delta f_{\lambda_1}', '\Delta f_{\lambda_2}','\Delta f_{\lambda_4}','\Delta f_{\lambda_5}'},'FontSize',16)
ylim([-0.1 0.1])
set(gcf,'position',[1327         841         673         497])
saveas(gcf,'figures/lambda_shift_stability_C','epsc')
saveas(gcf,'figures/lambda_shift_stability_C','png')




dX1 = X4 - X2;
dX2 = X4 - X3;
dX3 = X4 - X5;
dX4 = X4 - X6;

figure;
plot(T2, dX1,'LineWidth',1.5)
hold on
plot(T2, dX2,'LineWidth',1.5)
plot(T2, dX3,'LineWidth',1.5)
plot(T2, dX4,'LineWidth',1.5)
legend({'\Delta f_{SpanReg,\lambda_1}', '\Delta f_{SpanReg,\lambda_2}','\Delta f_{SpanReg,\lambda_4}','\Delta f_{SpanReg,\lambda_5}'},'FontSize',16)
ylim([-0.1 0.1])
set(gcf,'position',[1327         841         673         497])
saveas(gcf,'figures/lambda_shift_stability_D','epsc')
saveas(gcf,'figures/lambda_shift_stability_D','png')


% set(gcf,'position',[88          43        1698         937])
% saveas(gcf,'figures/lambda_shift_stability','epsc')
% saveas(gcf,'figures/lambda_shift_stability','png')
%%
% compute the norms
[norm(dx1) norm(dx2) norm(dx3) norm(dx4)]
[norm(dX1) norm(dX2) norm(dX3) norm(dX4)]