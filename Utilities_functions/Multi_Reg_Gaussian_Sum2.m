function [f_rec,alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum2(dat_noisy,id_OPT, adj_lambda, Gaus_info,Lambda)
%% this function uses the idea of MultiReg based on recovered lambda from GCV, which
% is known to have good recoverie.
% Lambda = Gaus_info.Lambda;
nLambda = length(Lambda);
nGaus = size(Gaus_info.LGBs,2);
A = Gaus_info.A;
m = size(A,2);
n = size(A,1);
f_unknown_L2 = zeros(m,nLambda);
T2 = Gaus_info.T2;
TE = Gaus_info.TE;
rho_L2 = zeros(nLambda,1);
eta_L2 = zeros(nLambda,1);SNR = 500;
n_run = 20;

%% A blend of offline computation
LGBs = Gaus_info.LGBs;
nc = size(LGBs,2);
if id_OPT == nLambda
    id_MultiReg = (id_OPT - adj_lambda):id_OPT;
elseif id_OPT == 1
    id_MultiReg = id_OPT:(id_OPT + adj_lambda);
else
    id_MultiReg = (id_OPT - adj_lambda):(id_OPT + adj_lambda);
end

lambda_MultiReg = Lambda(id_MultiReg);
%%
nLambda_temp = length(lambda_MultiReg);
beta_L2 = zeros(nc ,nLambda_temp);

L2_store = zeros(nc ,m,nLambda_temp);
options = optimoptions('quadprog','Display','off');

err_gaus_l2 = zeros(n_run,nc );
for k = 1:n_run
    
    new_L2_store = zeros(nc ,m,nLambda_temp);
    new_beta_store = zeros(nc, nLambda_temp);
    
%     new_ETA_L2 = zeros(nc ,nLambda);
    parfor i = 1:nc
        % calculate true exp.decaying signals
        dat_noiseless = A*LGBs(:,i);
        % add noise
        dat_noisy_ = dat_noiseless + max(abs(dat_noiseless))/SNR*randn(n,1);
        
        X_L2 = zeros(m,nLambda_temp);
        rho_temp = zeros(nLambda_temp,1);
        eta_temp = zeros(nLambda_temp,1);
        % for each noisy signals, calculate its L2 regularized solutions
        for kk = 1:length(lambda_MultiReg)
            [X_L2(:,kk),rho_temp(kk),eta_temp(kk)] = nonnegtik_hnorm(A,dat_noisy_,lambda_MultiReg(kk),'0');
        end
        
        % store the solutions
        new_L2_store(i,:,:) = X_L2;
        [beta_temp,err_gaus_l2(k,i)] = quadprog(X_L2'*X_L2, -X_L2'*LGBs(:,i),[],[],[],[],zeros(nLambda_temp,1),[],[],options);
        err_gaus_l2(k,i) = err_gaus_l2(k,i) + 1/2*LGBs(:,i)'*LGBs(:,i);
        new_beta_store(i,:) = beta_temp';
    end
    L2_store = L2_store + new_L2_store;
    beta_L2 = beta_L2 + new_beta_store;
    
end

L2_store = L2_store/n_run;
beta_L2 = beta_L2/n_run;
err_gaus_l2 = mean(err_gaus_l2);
%%
% store the averaged values
Gaus_info_temp.L2_store = L2_store;
Gaus_info_temp.beta_L2 = beta_L2;
Gaus_info_temp.LGBs = LGBs;
Gaus_info_temp.err_gaus_L2 = err_gaus_l2;
Gaus_info_temp.Lambda = lambda_MultiReg;
Gaus_info_temp.A = A;
Gaus_info_temp.T2 = T2;
Gaus_info_temp.TE = TE;
Gaus_info_temp.SNR = SNR;
%%
[f_rec,alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum1(dat_noisy,Gaus_info_temp);

end
