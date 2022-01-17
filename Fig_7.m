%% 
% This script generates brain images in terms of differences between ground
% truth
% 
% Chuan Bi
clc;clear;close all;
%% Load figure data
NESMA_results = load('brain_data/NESMA_filtered/sol_struct_brain_f49.mat').brain;
unfiltered_results = load('brain_data/unfiltered/sol_struct_brain_f49.mat').brain;
% load('f49/mask.mat');
% [p,q,r,s]=size(data);
load('brain_data/mask_2.mat')

SAVE = 1;
T2 = unfiltered_results.Gaus_info.T2;
Myelin_idx = find(T2<=40);

%% 

filtered_DP = NESMA_results.MWF_DP;
filtered_MR = NESMA_results.MWF_MR;
filtered_LS = NESMA_results.MWF_LS;
filtered_DP = BW.*imrotate(filtered_DP,275);
filtered_MR = BW.*imrotate(filtered_MR,275);
filtered_LS = BW.*imrotate(filtered_LS,275);
unfiltered_DP = unfiltered_results.MWF_DP;
unfiltered_MR = unfiltered_results.MWF_MR;
unfiltered_LS = unfiltered_results.MWF_LS;
unfiltered_DP = BW.*imrotate(unfiltered_DP,275);
unfiltered_MR = BW.*imrotate(unfiltered_MR,275);
unfiltered_LS = BW.*imrotate(unfiltered_LS,275);

% remove skulls and other unimportant regions
fig_filtered_DP = filtered_DP(60:250,80:240);
fig_filtered_MR = filtered_MR(60:250,80:240);
fig_filtered_LS = filtered_LS(60:250,80:240);
fig_unfiltered_DP = unfiltered_DP(60:250,80:240);
fig_unfiltered_MR = unfiltered_MR(60:250,80:240);
fig_unfiltered_LS = unfiltered_LS(60:250,80:240);


%% plot individual figures
figure;
ax(1) = subplot(2,3,1);
imagesc(fig_filtered_MR,[0,0.2]);axis off;axis equal;
title('SpanReg: high SNR', 'FontSize', 20)
ax(2) = subplot(2,3,2);
imagesc(fig_filtered_DP,[0,0.2]);axis off;axis equal;
title('DP: high SNR', 'FontSize', 20)
ax(3) = subplot(2,3,3);
imagesc(fig_filtered_LS,[0,0.2]);axis off;axis equal;
title('NNLS: high SNR', 'FontSize', 20)
ax(4) = subplot(2,3,4);
imagesc(fig_unfiltered_MR,[0,0.2]);axis off;axis equal;
title('SpanReg: low SNR', 'FontSize', 20)
ax(5) = subplot(2,3,5);
imagesc(fig_unfiltered_DP,[0,0.2]);axis off;axis equal;
title('DP: low SNR', 'FontSize', 20)
ax(6) = subplot(2,3,6);
imagesc(fig_unfiltered_LS,[0,0.2]);axis off;axis equal;
title('NNLS: low SNR', 'FontSize', 20)
h=colorbar;
% [left, bottom, width, height]
set(h, 'Location', 'south','Position',[0.1 0.05 0.84 0.05],'TickLabels',...
    {'0%','2%','4%','6%','8%','10%','12%','14%','16%','18%','20%'});
set(gcf,'Position',[440           5        1049         793]);
Name = 'figures/recoveries_comparison';
if SAVE == 1
    saveas(gcf,Name,'epsc')
    saveas(gcf,Name,'png')
end
%% calculate difference
diff_DP = fig_filtered_DP - fig_unfiltered_DP;
diff_MR = fig_filtered_MR - fig_unfiltered_MR;
diff_LS = fig_filtered_LS - fig_unfiltered_LS;

disp_DP = sprintf('The scaled error between low and high SNR for DP is %0.5f',sum(sum(abs(diff_DP)))/sum(sum(abs(fig_filtered_DP))));
disp_MR = sprintf('The scaled error between low and high SNR for MR is %0.5f',sum(sum(abs(diff_MR)))/sum(sum(abs(fig_filtered_MR))));
disp_LS = sprintf('The scaled error betweenlow and high SNR for NNLS is %0.5f',sum(sum(abs(diff_LS)))/sum(sum(abs(fig_filtered_LS))));
disp(disp_DP)
disp(disp_MR)
disp(disp_LS)

