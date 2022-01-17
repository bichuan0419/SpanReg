function Gaus_info = generate_gaussian_regs_L2_old(A,T2,TE,SNR,n_run,reg_param_lb,reg_param_ub,N_reg,Nc,cmin,cmax,...
    sigma_min,sigma_max)

% Given the true model Ax = b, compute the L2_regularized solutions by running
% multiple times of regularizations with noise added to each Gaussian basis
% functions.

% Input:
% A: Discrete Laplace Transform Matirx
% SNR: defines noise level by max(b)/SNR, additive white noise
% n_run: number of running times to be averaged to get mean values of
%        the outputs.
% reg_param_lb: lower bound for regularization parameters
% reg_param_ub: upper bound for regularization parameters
% nc: number of means in the Gaussian basis
% nsigma: number of standard deviations for the Gaussian basis
% cmin: lower bound for means for Gaussian basis
% cmax: upper bound for means for Gaussian basis
% sigma_min: lower bound for stds for Gaussian basis
% sigma_max: upper bound for stds for Gaussian basis
% T2: discrete T2 values

% Output: Matlab Structure named "Gaus_info", which contains
% beta: the beta values(coefficients) for each Gaussian basis
%           represented by their regularized solutions after averaging.
% L2_store: averaged L2_regularized solutions for all Gaussian basis
%           functions.
% LGBs: Gaussian basis.
% RHO_L2: averaged model error for each regularization for each Gaussian
%         basis.
% err_gaus_l2: averaged misfit for combination of regularized solutions to
%              represent the true Gaussian basis.

% Created by Chuan Bi, 03/22/2019
%% setting up Gaussian basis for computations
% Gaussian dictionary set up
LGBs = Gaussian_basis(T2,cmin,cmax,Nc,sigma_min,sigma_max);
options = optimoptions('quadprog','Display','off');
% reshape Gaussian dictionaries so that rows stands for Gaussian basis and
% columns stands for T2 values.
m = length(T2);
n = size(A,1);

% setting up regularization parameters
Lambda = logspace(reg_param_lb,reg_param_ub,N_reg);

nLambda = length(Lambda);
nc = sum(Nc);
%% Offline computation starts
% initializations
beta_L2 = zeros(nc ,nLambda);

L2_store = zeros(nc ,m,nLambda);

err_gaus_l2 = zeros(n_run,nc );
for k = 1:n_run
    
    new_L2_store = zeros(nc ,m,nLambda);
    new_beta_store = zeros(nc, nLambda);
    
%     new_ETA_L2 = zeros(nc ,nLambda);
    parfor i = 1:nc
        % calculate true exp.decaying signals
        dat_noiseless = A*LGBs(:,i);
        % add noise
        dat_noisy_ = dat_noiseless + max(dat_noiseless)/SNR*randn(n,1);
        
        X_L2 = zeros(m,nLambda);
        rho_temp = zeros(nLambda,1);
        eta_temp = zeros(nLambda,1);
        % for each noisy signals, calculate its L2 regularized solutions
        for kk = 1:length(Lambda)
            [X_L2(:,kk),rho_temp(kk),eta_temp(kk)] = nonnegtik_hnorm(A,dat_noisy_,Lambda(kk),'0');
        end
        
        % store the solutions
        new_L2_store(i,:,:) = X_L2;
        [beta_temp,err_gaus_l2(k,i)] = quadprog(X_L2'*X_L2, -X_L2'*LGBs(:,i),[],[],[],[],zeros(nLambda,1),[],[],options);
        err_gaus_l2(k,i) = err_gaus_l2(k,i) + 1/2*LGBs(:,i)'*LGBs(:,i);
        new_beta_store(i,:) = beta_temp';
    end
    L2_store = L2_store + new_L2_store;
    beta_L2 = beta_L2 + new_beta_store;
    
end

L2_store = L2_store/n_run;
beta_L2 = beta_L2/n_run;
err_gaus_l2 = mean(err_gaus_l2);



% store the averaged values
Gaus_info.L2_store = L2_store;
Gaus_info.beta_L2 = beta_L2;
Gaus_info.LGBs = LGBs;
Gaus_info.err_gaus_L2 = err_gaus_l2;
Gaus_info.Lambda = Lambda;
Gaus_info.A = A;
Gaus_info.T2 = T2;
Gaus_info.TE = TE;
Gaus_info.SNR = SNR;














end