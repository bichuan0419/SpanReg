%%
% this script is used to generate figures from the simulated human brain
% two SNR levels are considered: high: 800 and low: 200
%
% Chuan Bi
clc;clear;close all

%% load solved data
high_SNR_brain = load('brain_data/high_SNR/sol_struct_brain_f49.mat').brain;
low_SNR_brain = load('brain_data/low_SNR/sol_struct_brain_f49.mat').brain;
GT_brain = load('brain_data/ground_truth.mat').GT;
% load('f49/mask.mat');
% [p,q,r,s]=size(data);
load('brain_data/mask_2.mat')
%% get parameters;
SAVE = 1;

TE = high_SNR_brain.Gaus_info.TE;
T2 = high_SNR_brain.Gaus_info.T2;
A = high_SNR_brain.Gaus_info.A;
Myelin_idx = find(T2<=40);

%% get MWFs
high_MWF_DP = high_SNR_brain.MWF_DP;
high_MWF_LS = high_SNR_brain.MWF_LS;
high_MWF_MR = high_SNR_brain.MWF_MR;
low_MWF_DP = low_SNR_brain.MWF_DP;
low_MWF_LS = low_SNR_brain.MWF_LS;
low_MWF_MR = low_SNR_brain.MWF_MR;

% ground truth
[p,q,~] = size(GT_brain);
MWF_GT = zeros(p,q);

for j = 1:p
    for k = 1:q
        %GT
        temp = GT_brain(j,k,:);
        cur_T2_GT = temp(:);
        WF_GT = cumsum(cur_T2_GT);
        MWF_GT(j,k) = WF_GT(Myelin_idx(end));
    end
end
%% plot GT
fig_GT = BW.*imrotate(MWF_GT,275);fig_GT = fig_GT(60:250,80:240);
figure;
imagesc(fig_GT,[0,0.2]);axis off;axis equal;
title('Reference WMF map', 'FontSize', 27,'FontWeight','bold')
% h=colorbar;
% % [left, bottom, width, height]
% set(h,'TickLabels',...
%     {'0%','2%','4%','6%','8%','10%','12%','14%','16%','18%','20%'});
Name = 'figures/sim_GT';
if SAVE == 1
    saveas(gcf,Name,'epsc')
    saveas(gcf,Name,'png')
end
%% plot recoveries
fig_DP_high = BW.*imrotate(high_MWF_DP,275);fig_DP_high = fig_DP_high(60:250,80:240);
fig_MR_high = BW.*imrotate(high_MWF_MR,275);fig_MR_high = fig_MR_high(60:250,80:240);
fig_LS_high = BW.*imrotate(high_MWF_LS,275);fig_LS_high = fig_LS_high(60:250,80:240);

fig_DP_low = BW.*imrotate(low_MWF_DP,275);fig_DP_low = fig_DP_low(60:250,80:240);
fig_MR_low = BW.*imrotate(low_MWF_MR,275);fig_MR_low = fig_MR_low(60:250,80:240);
fig_LS_low = BW.*imrotate(low_MWF_LS,275);fig_LS_low = fig_LS_low(60:250,80:240);

figure;
ax(1) = subplot(2,3,1);
imagesc(fig_MR_high,[0,0.2]);axis off;axis equal;
title('SpanReg: high SNR', 'FontSize', 20)
ax(2) = subplot(2,3,2);
imagesc(fig_DP_high,[0,0.2]);axis off;axis equal;
title('DP: high SNR', 'FontSize', 20)
ax(3) = subplot(2,3,3);
imagesc(fig_LS_high,[0,0.2]);axis off;axis equal;
title('NNLS: high SNR', 'FontSize', 20)
ax(4) = subplot(2,3,4);
imagesc(fig_MR_low,[0,0.2]);axis off;axis equal;
title('SpanReg: low SNR', 'FontSize', 20)
ax(5) = subplot(2,3,5);
imagesc(fig_DP_low,[0,0.2]);axis off;axis equal;
title('DP: low SNR', 'FontSize', 20)
ax(6) = subplot(2,3,6);
imagesc(fig_LS_low,[0,0.2]);axis off;axis equal;
title('NNLS: low SNR', 'FontSize', 20)
h=colorbar;
% [left, bottom, width, height]
set(h, 'Location', 'south','Position',[0.1 0.05 0.84 0.05],'TickLabels',...
    {'0%','2%','4%','6%','8%','10%','12%','14%','16%','18%','20%'});
set(gcf,'Position',[440           5        1049         793]);
Name = 'figures/recoveries';
if SAVE == 1
    saveas(gcf,Name,'epsc')
    saveas(gcf,Name,'png')
end
%% calculate difference

diff_DP_high = fig_DP_high - fig_GT;
diff_MR_high = fig_MR_high - fig_GT;
diff_LS_high = fig_LS_high - fig_GT;

diff_DP_low = fig_DP_low - fig_GT;
diff_MR_low = fig_MR_low - fig_GT;
diff_LS_low = fig_LS_low - fig_GT;

disp_DP_high = sprintf('The scaled error between DP and GT at high SNR is %0.5f',sum(sum(abs(diff_DP_high)))/sum(sum(abs(fig_GT))));
disp_MR_high = sprintf('The scaled error between MR and GT at high SNR is %0.5f',sum(sum(abs(diff_MR_high)))/sum(sum(abs(fig_GT))));
disp_LS_high = sprintf('The scaled error between LS and GT at high SNR is %0.5f',sum(sum(abs(diff_LS_high)))/sum(sum(abs(fig_GT))));
disp_DP_low = sprintf('The scaled error between DP and GT at low SNR is %0.5f',sum(sum(abs(diff_DP_low)))/sum(sum(abs(fig_GT))));
disp_MR_low = sprintf('The scaled error between MR and GT at low SNR is %0.5f',sum(sum(abs(diff_MR_low)))/sum(sum(abs(fig_GT))));
disp_LS_low = sprintf('The scaled error between LS and GT at low SNR is %0.5f',sum(sum(abs(diff_LS_low)))/sum(sum(abs(fig_GT))));
disp(disp_DP_high)
disp(disp_MR_high)
disp(disp_LS_high)
disp(disp_DP_low)
disp(disp_MR_low)
disp(disp_LS_low)
