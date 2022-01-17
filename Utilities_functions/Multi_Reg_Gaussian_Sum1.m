function [f_rec_final,alpha_L2,F_info,C_L2] = Multi_Reg_Gaussian_Sum1(dat_noisy,Gaus_info)
%% Load information
Lambda = Gaus_info.Lambda;
nLambda = length(Lambda);
nGaus = size(Gaus_info.LGBs,2);
A = Gaus_info.A;
m = size(A,2);
n = size(A,1);
f_unknown_L2 = zeros(m,nLambda);
rho_L2 = zeros(nLambda,1);
eta_L2 = zeros(nLambda,1);

%% Online computation
for kk = 1:nLambda
    [f_unknown_L2(:,kk),rho_L2(kk),eta_L2(kk)] = nonnegtik_hnorm(A,dat_noisy,Lambda(kk),'0');
end

x_L2_weight = zeros(nGaus,nLambda);
    
Err_f_unknown_L2 = zeros(m,nLambda);

for j = 1:nLambda
    F_L2 = reshape(Gaus_info.L2_store(:,:,j),[],m)';
    Y_L2 = f_unknown_L2(:,j);
    %         obj_fit = @(xx) (F_L2*xx - Y_L2)'*(F_L2*xx - Y_L2);
    %         xx = fmincon(obj_fit,x0,[],[],[],[],lb,[],[],options);
    xx = lsqnonneg(F_L2,Y_L2);
    Err_f_unknown_L2(:,j) = F_L2*xx - Y_L2;
    x_L2_weight(:,j) = xx;
end

L2_store = Gaus_info.L2_store;
beta_L2 = Gaus_info.beta_L2;


L2store_reshaped_alpha = [];
L2store_reshaped_c = [];



for i = 1:nLambda
   L2store_reshaped_alpha = [L2store_reshaped_alpha,L2_store(:,:,i)'];
end
for i = 1:nGaus
    a = L2_store(i,:,:);
    L2store_reshaped_c = [L2store_reshaped_c, reshape(a,m,nLambda)];
end
% new_xweight has the form of [x_{11},x_{21},...,x_{M1},1|...|x_{1N},x_{2N},...,x_{MN},1]

new_xweight_L2 = x_L2_weight(:);


% new_beta has the form of [beta_{11},beta_{21},...,beta_{M1}|...|beta_{1N},beta_{2N},...,beta_{MN}]
beta_L2 = beta_L2';
new_beta_L2 = beta_L2(:);


% mat_alpha: each column is the contribution of g_{ji} from G_j to f_i
L2_mat_alpha = L2store_reshaped_alpha*diag(new_xweight_L2);


% mat_c: each column is the contribution of g_{ji} from 
L2_mat_c = L2store_reshaped_c*diag(new_beta_L2);

% % express in nLambda-1 variables
kron_alpha = kron(eye(nLambda),ones(nGaus,1));
% 
% % the way to get alpha in terms from unknown x = [alpha1...alphaN-1,c1...cM]
Ind_alpha = [eye(nLambda), zeros(nLambda,nGaus)];

% Now to get coefficients c
kron_c = kron(eye(nGaus),ones(nLambda,1));

% the way to get c in terms from unknown x = [alpha1...alphaN-1,c1...cM]
Ind_c = [zeros(nGaus,nLambda), eye(nGaus)];


% Now form linear ls problem
Aleft = L2_mat_alpha*kron_alpha*Ind_alpha - L2_mat_c*kron_c*Ind_c;
bright = zeros(m,1);
% precompute int_g_{i} for each g_i
% G = A(1,:)*Gaus_info.LGBs;
Aeq = [zeros(1,nLambda) ones(1,nGaus)];
beq = 1;
lb = zeros(nLambda + nGaus,1);
% options = optimoptions('quadprog','OptimalityTolerance',1e-16,'ConstraintTolerance',1e-16,'Display','off');
options = optimoptions('quadprog','Display','off');
H = Aleft'*Aleft;
% A_ineq = [ones(1,nLambda) zeros(1,nGaus)];
% b_ineq = 1;
s = quadprog(H,zeros(nLambda+nGaus,1),[],[],Aeq,beq,lb,[],[],options);
C_L2 = s(nLambda+1:end);
alpha_L2 = s(1:nLambda);

%% find the best representation
% res_alpha = norm(dat_noisy - A*f_unknown_L2*alpha_L2);
% res_c = norm(dat_noisy - G*C_L2);

f_rec_final = f_unknown_L2*alpha_L2;

%% Getting results and save data
f_rec = f_unknown_L2*alpha_L2;
F_info.alpha_L2 = alpha_L2;
F_info.f_unknown_L2 = f_unknown_L2;
F_info.rho_L2 = rho_L2;
F_info.eta_L2 = eta_L2;
F_info.f_rec = f_rec;
F_info.C_L2 = C_L2;


end
