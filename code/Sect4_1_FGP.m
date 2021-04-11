%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains an example to implement FGP, kriging with Matern 
% covariance function to obtain the results in Section 4.1 of 
% manuscript.
%
% The basis functions here are chosen using the R package FRK with version 0.2.1. 
% They are saved in the file named basis_info_2Res_2500grids.txt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% simulate from exponential covariance function 
clear; clc;

addpath('./FGP')

N = 2500;

nx = 50;
ny = 50;
s1 = 0;
s2 = 50;
L0 = s2 - s1;
xg = linspace(s1, s2, nx)';
yg = linspace(s1, s2, ny)';
[x, y] = meshgrid(xg, yg);
x = reshape(x, nx*ny, 1);
y = reshape(y, nx*ny, 1);
total_locs = [x, y];

dists = pdist2(total_locs,total_locs);
sig2eps = 4;
sig2=16;
phi=3;
mat_cov = sig2*exp(-dists/phi);

for i=1:N
	mat_cov(i,i) = sig2;
end

%L = chol(mat_cov+sig2eps*speye(N), 'lower');
L = chol(mat_cov, 'lower');

%%% specify locations for estimation and prediction
% 10% random selected location for out-of-sample prediction
% 90% as observation
MBD = [];
Remain = setdiff(1:N, MBD);

rng(123);	
% ind = unidrnd(length(Remain), ceil(length(Remain)*0.1), 1); 
ind = datasample(Remain, floor(length(Remain)*0.1), 'Replace', false); 

MAR = ind;	
MAR = MAR(:);
OBS = setdiff(Remain, MAR);
OBS = OBS(:);
locs = total_locs(OBS,:);
% pred_locs = total_locs;
% predID = 1:N';
pred_locs = total_locs(MAR,:);
predID = MAR;



basis_info=load('basis_info_2Res_2500grids.txt');
grid_all1 = basis_info(:, 1:2);
threshold1 = basis_info(:, 3);
nS1 = ones(length(threshold1), 1);



%% construct proximity matrix
dtemp = pdist2(total_locs, total_locs);
[indi, indj] = find(dtemp);
d0 = 1.01*sqrt((total_locs(1,1)-total_locs(2,1))^2 + (total_locs(1,2)-total_locs(2,2))^2); % first order neighbors
[indi, indj] = find(dtemp>0 & dtemp<=d0);
H = sparse(indi, indj, ones(length(indi),1), length(x), length(x));
summary(dataset(sum(H,2)))
reord = symrcm(H);
opts.tol=1e-5;
large = eigs(H(reord,reord), 1, 'LA', opts);
opts.tol=1e-3;
small = eigs(H(reord,reord), 1, 'SA', opts);

Hall.H = H;
Hall.large = large;
Hall.small = small;

% create indicator matrix
Ao = sparse(size(locs,1), length(x));
for k = 1:size(total_locs,1)
	Ao((locs(:,1)==total_locs(k,1)) & ...
		locs(:,2)==total_locs(k,2), k) = 1;
end
A{1}=Ao;

m = size(pred_locs, 1);
Ap = sparse(m, size(total_locs,1));
for k = 1:size(total_locs,1)
	Ap((pred_locs(:,1)==total_locs(k,1)) & ...
    pred_locs(:,2)==total_locs(k,2), k) = 1;
end
A{2}=Ap;

%%% pre-allocation
nsim = 15;
ind = cell(nsim, 1);
Y = cell(nsim, 1);
Y_test = cell(nsim, 1);
Z_test = cell(nsim, 1);


MSPE1_MAR = zeros(nsim, 1);
MSPE2_MAR = zeros(nsim, 1);
MSPE3_MAR = zeros(nsim, 1);


crps1_MAR = zeros(nsim, 1);
crps2_MAR = zeros(nsim, 1);
crps3_MAR = zeros(nsim, 1);




tic;
for i = 1:nsim


%i=1;
Ysim = L*randn(N,1);
Zsim = Ysim + sqrt(sig2eps)*randn(N,1);
Y{i} = Zsim(OBS);
Y_test{i} = Ysim(predID);
%Z_test{i} = Ysim(predID);



%%%%%%%%% Exact Kriging %%%%%%
ts1=tic;
ind_locs =OBS;
cross_cov = mat_cov(ind_locs, predID);
cov_obs = mat_cov(ind_locs, ind_locs) + sig2eps*speye(length(ind_locs));
Lcov = chol(cov_obs, 'lower');
LcovY = Lcov\Y{i};
Lcros = Lcov\cross_cov;
pred1 = Lcros'*LcovY;

cov_pred = mat_cov(predID, predID);
sig2_pred1 = diag(cov_pred) - diag(Lcros'*Lcros);

% MSPE1_MAR(i) = sum((Y_test{i}(MAR) - pred1{i}(MAR)).^2)/length(MAR);
% crps1_MAR(i) = mean(CRPS(Y_test{i}(MAR), pred1{i}(MAR), sig2_pred1{i}(MAR).^0.5));

MSPE1_MAR(i) = mean((Y_test{i} - pred1).^2);
crps1_MAR(i) = mean(CRPS(Y_test{i}, pred1, sig2_pred1.^0.5));

te1(i)=toc(ts1);



%%%%%%%%% Kriging with fitted statioanry Matern (exp covariance) %%%%%%%%%%

%% estimation
Ytilde = Y{i};
[rho, sig2_variance] = MLE_Matern(Ytilde, locs, 0.5, sig2eps);
[pred_notrend, sig2_pred2] = spKrigingMatern(Ytilde, locs, ...
	pred_locs, 0.5, rho, sig2_variance, sig2eps);
pred2 = pred_notrend;

MSPE2_MAR(i) = mean((Y_test{i} - pred2).^2);
crps2_MAR(i) = mean(CRPS(Y_test{i}, pred2, sig2_pred2.^0.5));



%%%%%%%%%%%%%%%%%% model 2: FGP %%%%%%%%%%%%%%%%%%%%%%%

[pred3, sig2_pred3, tau2, gamma] = wrapFGP(locs, ...
	Y{i}, [], pred_locs, [], grid_all1, threshold1, Hall, A, sig2eps);

MSPE5_MAR(i) = mean((Y_test{i} - pred3).^2);
crps5_MAR(i) = mean(CRPS(Y_test{i}, pred3, sig2_pred3.^0.5));


end
timing = toc;


%% EK MK FGP
out1=[mean(MSPE1_MAR), mean(MSPE2_MAR), mean(MSPE3_MAR)];
out2=[std(MSPE1_MAR), std(MSPE2_MAR), std(MSPE3_MAR)];

MSPE=vertcat(out1,out2)


mu_CRPS_MAR=[mean(crps1_MAR), mean(crps2_MAR), mean(crps3_MAR)];
std_CRPS_MAR=[std(crps1_MAR), std(crps2_MAR), std(crps3_MAR)];

crps=vertcat(mu_CRPS_MAR,std_CRPS_MAR)




