function [avg_MDL_err, avg_MDL_err_DP] = heatmap_unequal_width_DP(Gaus_info, show, n_sim)
%% Generate a sequence of Simulation
T2 = Gaus_info.T2;
TE = Gaus_info.TE;
A = Gaus_info.A;
npeaks = 2;
% ratio of peak separation
rps = linspace(1,4,5)';
mps = rps/2;
nrps = length(rps);
T2_left = 30*ones(nrps,1);
T2_mid = T2_left.*mps;
T2_right = T2_left.*rps;
T2mu = [T2_left T2_right];
% T2mu = [T2_left T2_mid T2_right];
% Uniform Gaussian SD
nsigma = 5;
unif_sigma = linspace(2,5,nsigma)';
diff_sigma = [unif_sigma, 3*unif_sigma];
f_coef = ones(npeaks,1);
%%
LGBs = Gaus_info.LGBs;
Lambda = Gaus_info.Lambda;
Lambda_1 = logspace(log10(Lambda(1)), log10(Lambda(end)), length(Lambda));
% Lambda_1 = Lambda;
% if size(Lambda_1,1)<size(Lambda_1,2)
%     Lambda_1 = Lambda_1';
% end
nLambda_1 = length(Lambda_1);
nLambda = length(Lambda);
nGaus = size(Gaus_info.LGBs,2);
[n,m] = size(Gaus_info.A);
avg_MDL_err = zeros(nsigma,nrps);
avg_MDL_err_DP = zeros(nsigma,nrps);

error_multi_reg = zeros(nsigma,nrps,n_sim);
error_DP = zeros(nsigma,nrps,n_sim);


SNR = Gaus_info.SNR;
%% create dataset to save
MultiReg_data = zeros(n_sim,nsigma,nrps, m);
DP_data = zeros(n_sim,nsigma,nrps, m);
IdealModel_data = zeros(n_sim,nsigma,nrps, m);
%%
for l = 1:n_sim
    MDL_err = zeros(nsigma,nrps);
    MDL_err_dp = zeros(nsigma,nrps);

    for i = 1:nsigma
        sigma_i = diff_sigma(i,:);
%         sigma_i1 = unif_sigma(nsigma-i);

        if show == 0
            parfor j = 1:length(rps)
                % setting up the distributions
                p = zeros(npeaks,m);
                T2mu_sim = T2mu(j,:);
                for ii = 1:npeaks
                    p(ii,:) = normpdf(T2,T2mu_sim(ii),sigma_i(ii));
                end
    %             IdealModel_weighted = p'*(f_coef.*rand(npeaks,1));
                IdealModel_weighted = p'*(f_coef)/npeaks;
                
                dat_noiseless = A*IdealModel_weighted;
                dat_noisy = dat_noiseless + max(abs(dat_noiseless))/SNR*randn(length(TE),1);


    %             figure;
    %             plot(TE,dat_noisy,'LineWidth',2);
    %             ylim([-0.1,1.2])
    %             xlabel('TE','FontSize',16)
    %             ylabel('y_{ob}','FontSize',16)
    %             title('Synthetic Data','FontSize',16)
    %             
                disp(['evaluating ',num2str(i),'-th sigma and ',num2str(j),'-th ratio peak separation test'])
                % online computation
                %% DP
                [f_rec_dp,lambda_dp] = discrep_L2(dat_noisy,A,SNR,Lambda);
                %% online computation
                [f_rec,alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum1(dat_noisy,Gaus_info);
                %% save to dataset
                MultiReg_data(l,i,j,:) = f_rec;
                DP_data(l,i,j,:) = f_rec_dp;
                
                %%
                true_norm = norm(IdealModel_weighted);
                MDL_err(i,j) = norm(IdealModel_weighted - f_rec)/true_norm;
                MDL_err_dp(i,j) = norm(IdealModel_weighted - f_rec_dp)/true_norm;

                error_multi_reg(i,j,l) = MDL_err(i,j);
                error_DP(i,j,l) = MDL_err_dp(i,j);
                

            end
        else
            for j = 1:length(rps)
                % setting up the distributions
                p = zeros(npeaks,m);
                T2mu_sim = T2mu(j,:);
                for ii = 1:npeaks
                    p(ii,:) = normpdf(T2,T2mu_sim(ii),sigma_i(ii));
                end
    %             IdealModel_weighted = p'*(f_coef.*rand(npeaks,1));
                IdealModel_weighted = p'*(f_coef)/npeaks;
                dat_noiseless = A*IdealModel_weighted;
                dat_noisy = dat_noiseless + max(abs(dat_noiseless))/SNR*randn(length(TE),1);


    %             figure;
    %             plot(TE,dat_noisy,'LineWidth',2);
    %             ylim([-0.1,1.2])
    %             xlabel('TE','FontSize',16)
    %             ylabel('y_{ob}','FontSize',16)
    %             title('Synthetic Data','FontSize',16)
    %             
                disp(['evaluating ',num2str(i),'-th sigma and ',num2str(j),'-th ratio peak separation test'])
                % online computation
                %% DP
                [f_rec_dp,lambda_dp] = discrep_L2(dat_noisy,A,SNR,Lambda);
                %% online computation
                [f_rec,alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum1(dat_noisy,Gaus_info);
                %% save to dataset
                MultiReg_data(l,i,j,:) = f_rec;
                DP_data(l,i,j,:) = f_rec_dp;
                IdealModel_data(l,i,j,:) = IdealModel_weighted;
                %%
                true_norm = norm(IdealModel_weighted);
                MDL_err(i,j) = norm(IdealModel_weighted - f_rec)/true_norm;
                MDL_err_dp(i,j) = norm(IdealModel_weighted - f_rec_dp)/true_norm;

                error_multi_reg(i,j,l) = MDL_err(i,j);
                error_DP(i,j,l) = MDL_err_dp(i,j);
                %%
%                 figure;
%                 set(gcf,'position',[284         538        1206         420]);
%                 subplot(1,2,1)
%                 plot(T2,IdealModel_weighted,'LineWidth',3,'color',[0,0,0]);
%                 hold on
%                 plot(T2,f_rec,':','LineWidth',3,'color',[1, 0, 0]);
%                 plot(T2, Gaus_info.LGBs*C_L2,'LineWidth',3,'color',[1, 1, 0])
%                 plot(T2,f_rec_dp,'--','LineWidth',3,'color',[0, 0, 1]);
% %                 legend({'True Dist.','Multi-Reg','DP'},'FontSize',20,'location','best')
%                 legend({'True Dist.','Multi-Reg \alpha','MR c','DP'},'FontSize',20,'location','best')
%                 xlabel('T2 Relaxation Time','FontSize',20,'FontWeight','bold');
%                 ylabel('Intensity','FontSize',20,'FontWeight','bold');
%                 title({'Comparison of Recovered Distributions'},'FontSize',20,'FontWeight','bold')
%                 subplot(1,2,2)
%                 plot(TE,A*IdealModel_weighted,'LineWidth',3,'color',[0,0,0]);
%                 hold on
%                 plot(TE,A*f_rec,':','LineWidth',3,'color',[1, 0, 0]);
%                 plot(TE,A*f_rec_dp,'--','LineWidth',3,'color',[0, 0, 1]);
%                 legend({'True Dist.','Multi-Reg','DP'},'FontSize',20,'location','best')
%                 xlabel('TE','FontSize',20,'FontWeight','bold');
%                 ylabel('Amplitude','FontSize',20,'FontWeight','bold');
%                 title({'Comparison of Recovered Data'},'FontSize',20,'FontWeight','bold')
%                 drawnow
                
                
            end
        end
    end
    avg_MDL_err = avg_MDL_err + MDL_err;
    avg_MDL_err_DP = avg_MDL_err_DP + MDL_err_dp;

end
avg_MDL_err = avg_MDL_err/n_sim;
avg_MDL_err_DP = avg_MDL_err_DP/n_sim;

compare_DP = avg_MDL_err-avg_MDL_err_DP;

sol_strct_uneql.avg_MDL_err = avg_MDL_err;
sol_strct_uneql.avg_MDL_err_DP = avg_MDL_err_DP;
sol_strct_uneql.compare_DP = compare_DP;
sol_strct_uneql.Gaus_info = Gaus_info;
sol_strct_uneql.error_multi_reg = error_multi_reg;
sol_strct_uneql.error_DP = error_DP;

% save history
sol_strct_uneql.MultiReg_data = MultiReg_data;
sol_strct_uneql.DP_data = DP_data;
if show == 1
    sol_strct_uneql.IdealModel_data = IdealModel_data;
end

if show == 0
    FileName = 'sol_struct_uneql2_noshow_old_DP';
    matfile = fullfile('figure_data', FileName);
    save(matfile,'sol_strct_uneql')
else
    FileName = 'sol_struct_uneql2_show_old_DP';
    matfile = fullfile('figure_data', FileName);
    save(matfile,'sol_strct_uneql')
end


% figure;
% image(rps,unif_sigma,compare_DP,'CDataMapping','scaled')
% xlabel('Ratio of Peak Separation','Fontsize',18)
% ylabel('Gaussian \sigma','Fontsize',18)
% c = colorbar;
% c.Label.String = 'H_{MultiReg} - H_{DP}';
% c.Label.FontSize = 18;
% 
% title('Comparison with DP','Fontsize',18)
% caxis([-0.3 0.3])

end