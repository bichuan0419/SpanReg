clc;clear;close all
%% offline computation0.*
% rng(0)
run = 0;
load lambda_16_SNR_500_nrun_20_sigma_min_2_sigma_max_4_basis2_1604020lmbda_min-6lmbda_max1.mat
%%
show = 0;
n_sim = 20;
%%
[avg_MDL_err, avg_MDL_err_DP] = heatmap_unequal_width_DP(Gaus_info, show, n_sim);
difference = avg_MDL_err - avg_MDL_err_DP;
MultiReg_MDL = mean(mean((avg_MDL_err)));
DP_MDL = mean(mean((avg_MDL_err_DP)));
%% plot
load 'figure_data/sol_struct_uneql2_noshow_old_DP.mat'
rps = linspace(1,4,5)';
nsigma = 5;
unif_sigma = linspace(2,5,nsigma)';

compare_DP = sol_strct_uneql.compare_DP;
avg_MDL_err_DP = sol_strct_uneql.avg_MDL_err_DP;
avg_MDL_err = sol_strct_uneql.avg_MDL_err;


figure;
image(rps,unif_sigma,compare_DP,'CDataMapping','scaled')
% plot(0,0)
xlabel('Ratio of Peak Separation','Fontsize',42)
ylabel('Gaussian \sigma','Fontsize',42)
c = colorbar;
caxis([-.5 0.3])
c.FontSize = 24;
c.Location = 'southoutside';
annotation('textbox', [0.1, 0.02, 0.4, 0.05], 'string', 'SpanReg superior','Fontsize',36,'EdgeColor','none')
annotation('textbox', [0.65, 0.02, 0.95, 0.05], 'string', 'SpanReg inferior','Fontsize',36,'EdgeColor','none')
set(gcf,'position',[423          40        1260        1217])
title('SpanReg versus DP','Fontsize',55,'FontWeight','bold')
saveas(gcf,'figures/Heatmap_accuracy','epsc')
saveas(gcf,'figures/Heatmap_accuracy','png')
%%
figure;
image(rps,unif_sigma,avg_MDL_err,'CDataMapping','scaled')
xlabel('Ratio of Peak Separation','Fontsize',42)
ylabel('Gaussian \sigma','Fontsize',42)
c = colorbar;
caxis([0 1])
c.FontSize = 24;
c.Location = 'southoutside';
annotation('textbox', [0.1, 0.02, 0.4, 0.05], 'string', 'Better recovery','Fontsize',36,'EdgeColor','none')
annotation('textbox', [0.65, 0.02, 0.95, 0.05], 'string', 'Worse recovery','Fontsize',36,'EdgeColor','none')
set(gcf,'position',[423          40        1260        1217])
title('SpanReg','Fontsize',55,'FontWeight','bold')
saveas(gcf,'figures/Heatmap_MR','epsc')
saveas(gcf,'figures/Heatmap_MR','png')

%%
figure;
image(rps,unif_sigma,avg_MDL_err_DP,'CDataMapping','scaled')
xlabel('Ratio of Peak Separation','Fontsize',42)
ylabel('Gaussian \sigma','Fontsize',42)
c = colorbar;
caxis([0 1])
c.FontSize = 24;
c.Location = 'southoutside';
annotation('textbox', [0.1, 0.02, 0.4, 0.05], 'string', 'Better recovery','Fontsize',36,'EdgeColor','none')
annotation('textbox', [0.65, 0.02, 0.95, 0.05], 'string', 'Worse recovery','Fontsize',36,'EdgeColor','none')
set(gcf,'position',[423          40        1260        1217])
title('DP','Fontsize',55,'FontWeight','bold')
saveas(gcf,'figures/Heatmap_DP','epsc')
saveas(gcf,'figures/Heatmap_DP','png')

