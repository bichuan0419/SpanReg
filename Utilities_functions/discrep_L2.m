function [f_rec_dp, lambda] = discrep_L2(dat_noisy,A,SNR,Lambda)
% this code does conventional l2 regularization with positivity constraints
% on f based on the model
% y = argmin{ || y - Af ||_2^2 + lambda^2 || f ||_2^2 } such that f >= 0

%% Load information
nLambda = length(Lambda);
m = size(A,2);
n = size(A,1);

threshhold = 1.05*sqrt(n)*max(dat_noisy)/SNR;
lambda = Lambda(1);
diff_rho = -1;
X = [];
while diff_rho <= 0
    lambda = lambda*1.5;
    [X,rho] = nonnegtik_hnorm(A,dat_noisy,lambda,'0');
    diff_rho = sqrt(rho) - threshhold;
end
f_rec_dp = X;
end