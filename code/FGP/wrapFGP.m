	
% Author: Pulong Ma
% Date:  June 04th, 2015
% Last Modified by: Pulong Ma
% Last Modified Date: October 27, 2017
% Last Modified time: 10:31:08


% Purpose: 


function [pred, sig2_pred, tau2, gamma, beta, sig2eps] = wrapFGP(locs, ...
		 Y, Xo, pred_locs, Xp, grid_all, threshold, Hall, A, sig2eps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%

lat = locs(:, 1);
lon = locs(:, 2);

pred_lat = pred_locs(:, 1);
pred_lon = pred_locs(:, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% detrend the data %%%%%%%%%%%%%%%%%%%%%

% create the matrix of predictors 
n = length(Y);
if isempty(Xo) || isempty(Xp)
	Xo = [ones(n, 1)];	% constant trend
	%Xo = [ones(n, 1), lat, lon];
	%Xo = [ones(n, 1), lat, lat.^2];	% quadratic trend
	Xp = [ones(length(pred_lat), 1)];	
	%Xp = [ones(length(pred_lat), 1), pred_lat, pred_lat.^2];
	%Xp = [ones(length(pred_lat), 1), pred_lat, pred_lon];	
end



% create wendland basis function
S = bisquare(locs, grid_all, threshold);

% estimate parameters

%[K,tau2,gamma,T,epsilon,time] = EM_SRECAR(S,A{1},Y_tilde,Hall,sig2eps);
%disp(strcat('estimation error:', ' ', num2str(epsilon)));

if  (~exist('sig2eps', 'var'))
	[K,tau2,gamma,beta,sig2eps,T,epsilon,time] = SEM_FGP(S,A{1},Y,Xo,Hall,[]);
	disp(strcat('estimation error:', ' ', num2str(epsilon)));
	disp(strcat('time for parameter estimation: ', num2str(time), 'seconds'))
else
	[K,tau2,gamma,beta,~,T,epsilon,time] = SEM_FGP(S,A{1},Y,Xo,Hall,sig2eps);
	disp(strcat('estimation error:', ' ', num2str(epsilon)));
	disp(strcat('time for parameter estimation: ', num2str(time), 'seconds'))
end 
%%%%%%%%%%%%%%%%%%%%%%%% obtain predictions %%%%%%%%%%%%%%%%%%%
Y_tilde = Y - Xo*beta;
data = [lat lon Y_tilde, Xo];

[pred_notrend, sig2_pred] = FGP(data, pred_locs, Xp, K, tau2, gamma, ...
				Hall.H, A, grid_all, threshold, sig2eps);

% add back trend

pred = pred_notrend + Xp*beta;



