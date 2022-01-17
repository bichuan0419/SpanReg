clc;clear;close all
%% offline computation
rng(11)
run = 0;
load lambda_16_SNR_500_nrun_20_sigma_min_2_sigma_max_4_basis2_1604020lmbda_min-6lmbda_max1.mat
%%
show = 1;
n_sim = 1;
%%
[avg_MDL_err, avg_MDL_err_DP] = heatmap_unequal_width_DP(Gaus_info, show, n_sim);
difference = avg_MDL_err - avg_MDL_err_DP;
MultiReg_MDL = mean(mean((avg_MDL_err)));
DP_MDL = mean(mean((avg_MDL_err_DP)));
%% plot
load('figure_data/sol_struct_uneql2_show_old_DP.mat');
Gaus_info = sol_strct_uneql.Gaus_info;
T2 = Gaus_info.T2;
[n_sim,nsigma,nrps, m] = size(sol_strct_uneql.MultiReg_data);
MultiReg_data = reshape(sol_strct_uneql.MultiReg_data(1,:,:,:),nsigma,nrps,[]);
DP_data = reshape(sol_strct_uneql.DP_data(1,:,:,:),nsigma,nrps,[]);
IdealModel_data = reshape(sol_strct_uneql.IdealModel_data(1,:,:,:),nsigma,nrps,[]);
% plot all figures;
fig = figure;
set(gcf,'position',[680         152        1089         826])
for i = 1:nsigma
    for j = 1:nrps
        disp([nsigma*(i-1)+j])


        subplot(nsigma,5,5*(i-1)+j)
        hold on
        plot(T2,reshape(IdealModel_data(i,j,:),m,1),'LineWidth',2,'color',[0,0,0]);
        plot(T2,reshape(MultiReg_data(i,j,:),m,1),'LineWidth',2,'color',[1,0,0]);
        plot(T2,reshape(DP_data(i,j,:),m,1),'LineWidth',2,'color',[0,0,1]);
        hold off
        legend({'True','SpanReg','DP'})
%             xlabel('T2','FontSize',20,'FontWeight','bold');
        drawnow
    end
end
saveas(gcf,'figures/Recoveries_Comparision','epsc')
saveas(gcf,'figures/Recoveries_Comparision','png')
% title('Recoveries Comparison')