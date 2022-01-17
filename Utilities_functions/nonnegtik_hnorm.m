function [x,mdl_err,xnorm] = nonnegtik_hnorm(A,b,lambda,nm,R)
% this code deal with 
% (1) with Tikhonov Regularization equiped with Hp norm, p>=1
% (2) with nonnegative constrains on x
% (3) norm: string, 0<=norm<=2 is the H_p norm
% function minimizes |Ax - b|^2 + alpha|x|_{H^p}^2 for x_i >= 0 for all i
% Chuan Bi
% example: [x,L] = nonnegtik_hnorm(A,b,lambda,T2,'2')

n = size(A,2);

L0 = eye(n);

L1 = -diag(ones(n,1)) + diag(ones(n-1,1),1);
L1(end,:) = [];

L2 = (6*diag(ones(n,1))  - 4*diag(ones(n-1,1),1) - 4*diag(ones(n-1,1),-1) + ...
    diag(ones(n-2,1),-2) + diag(ones(n-2,1),2));

Aug_b = [b; zeros(n,1)];
switch nargin
    case 3
        L = L0;
    case 4
        switch nm
            case '0'
                L = L0;
                
            case '1'
                L = L0 + L1;
                
            case '2'
                L = L0 + L1 + L2;
                
            case '11'
                L = L1;
                
            case '22'
                L = L2;
        end
    case 5
        L = R;
        
end

Aug_A = [A;lambda*L];
x=lsqnonneg(Aug_A,Aug_b);

mdl_err = norm(A*x - b)^2;
xnorm = norm(L*x)^2;
