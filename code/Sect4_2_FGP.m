%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code reproduces the results based on FGP in Section 4.2 of 
% manuscript.
%
% The basis functions here are chosen using the R package FRK with version 0.1.8. 
% They are saved in the file named basis_info_2Res_10grids.txt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

addpath('./FGP')



f = @(s) exp(-(s/750-1).^2) + exp(-0.8*(s/750+1).^2) - 0.05*sin(8*(s/750+0.1));

nx = 100;
ny = 100;
s1 = -1500;
s2 = 1500;
L0 = s2 - s1;
xg = linspace(s1, s2, nx)';
yg = linspace(s1, s2, ny)';
[x, y] = meshgrid(xg, yg);
x = reshape(x, nx*ny, 1);
y = reshape(y, nx*ny, 1);
total_locs = [x, y];


%% read random numbers generated from R 
F = -10*f(x).*f(y);

sig2eps = 0.1*var(F);



%%% prediction locations 
n = nx*ny;
m = floor(n*0.1);

%flagMR = zeros(n,1);

flagMR = total_locs(:,1)>s1+0.25*L0 & total_locs(:,1)<s1+0.5*L0 & ...
 total_locs(:,2)>s1+0.375*L0 & total_locs(:,2)<s1+0.625*L0;

total_ID = 1:n;
total_ID = total_ID(:); 
MBD_ID = total_ID(flagMR==1);
RemainID = total_ID(flagMR==0);


rng(230);
indrand = datasample(RemainID, m, 'Replace', false);
ind_locs = setdiff(RemainID, indrand);

% dlmwrite('MAR_ind_sim2D.txt', indrand, 'delimiter', '\t');


predID = [MBD_ID;indrand];
pred_locs = total_locs(predID, :);

locs = total_locs(ind_locs, :);



basis_info=load('basis_info_2Res_10kgrids.txt');
grid_all1 = basis_info(:, 1:2);
threshold1 = basis_info(:, 3);



%%% create proximity matrix
% first order proximity matrix
dtemp = pdist2(total_locs, total_locs);
[indi, indj] = find(dtemp);
% d0 = 1.2*min(L0/nx, L0/ny);
d0 = 1.01*sqrt((total_locs(1,1)-total_locs(2,1))^2 + (total_locs(1,2)-total_locs(2,2))^2); % first order neighbors
[indi, indj] = find(dtemp>0 & dtemp<=d0);
H = sparse(indi, indj, ones(length(indi),1), size(total_locs,1), size(total_locs,1));
reord = symrcm(H);
large = eigs(H(reord,reord), 1, 'LA');
small = eigs(H(reord,reord), 1, 'SA');
Hall.H = H;
Hall.large = large;
Hall.small = small;
% clear dtemp indi indj;

%% create indicator matrix
Ao = sparse(size(locs,1), size(total_locs,1));
for k = 1:size(total_locs,1)
	Ao((locs(:,1)==total_locs(k,1)) & locs(:,2)==total_locs(k,2), k) = 1;
end
A{1}=Ao;

%pred_locs = total_locs;

m = size(pred_locs, 1);
Ap = sparse(m, size(total_locs,1));
for k = 1:size(total_locs,1)
	Ap((pred_locs(:,1)==total_locs(k,1)) & ...
    pred_locs(:,2)==total_locs(k,2), k) = 1;
end
A{2}=Ap;





%%% pre-alocation
a = 15;
pred = cell(a,1);
sig2pred = cell(a,1);

tau2 = zeros(a,1);
gamma = zeros(a,1);
beta = cell(a,1); % FGP
sig2eps_em = zeros(a,1);

% FGP
MSPE = zeros(a,1);
crps = zeros(a,1);


Y_sim = cell(a, 1);
Y_test = cell(a, 1);
Ytrue = cell(a, 1);
Z_sim = cell(a, 1);
Z_test = cell(a, 1);






rng(123);


for i=1:a

Ysim = F;

Y_sim{i} = Ysim;

Z_sim{i} = Y_sim{i} + sqrt(sig2eps)*randn(nx*ny,1);


%%%%%%%%% create data %%%%%%%%%%
Z = Z_sim{i};
% create data for estimation
Y = Z(ind_locs);


Y_test{i} = F(predID);





%%%%%%%%%% FGP model %%%%%%%%%%%
ts1 = tic;
[pred{i}, sig2pred{i}, tau2(i), gamma(i), beta{i}, sig2eps_em(i)] = wrapFGP(locs, ...
Y, [], pred_locs, [], grid_all1, threshold1, Hall, A, sig2eps);
te1(i) = toc(ts1);

MSPE(i) = mean((Y_test{i} - pred{i}).^2);
crps(i) = mean(CRPS(Y_test{i}, pred{i}, (sig2pred{i}).^0.5));




end



%[mean(MSPE), std(MSPE); mean(crps), std(crps)]
[mean(sqrt(MSPE)); std(sqrt(MSPE)); mean(crps); std(crps)]




