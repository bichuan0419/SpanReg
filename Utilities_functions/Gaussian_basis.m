function LGBs = Gaussian_basis(x,cmin,cmax,Nc,sigma_min,sigma_max)

% create a set of gaussian basis functions, which is represented
% by a matrix, whose columns are of size nx (or n(T2)) and the rows are
% of size sum(Nc)
% each column is one Gaussian distribution
% c stands for centers
% sigma stands for standard deviations

%%
% 
nsigma = length(Nc);
sigma = linspace(sigma_min,sigma_max,nsigma); % set of standard deviations

c = cell(nsigma,1);

for k = 1:nsigma
    c{k} = linspace(cmin + 3*sigma(k),cmax - 3*sigma(k),Nc(k)); % set of centers for each sigma
end


[a,b] = size(x);
if a<b; x = x';end % x => Nx1 vector
% nc = sum(Nc);
LGBs = [];
for j = 1:nsigma
    for i = 1:length(c{j})
        gaus = normpdf(x,c{j}(i),sigma(j));
        LGBs = [LGBs, gaus];
    end
end