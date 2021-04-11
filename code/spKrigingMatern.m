	
% Author: Pulong Ma University of Cincinnati
% Date: September 30, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: October 05, 2015
% Last Modified time: 20:59:37


% Purpose: Simple Kriging with Matern covariance function



function [pred_notrend, sig2pred] = spKrigingMatern(z, locs, pred_locs, nu, rho, sig2, sig2eps)
%% Input Arguments:
% z: detailed residuals
% locs: locations of observations
% pred_locs: prediction locations
% nu: smoothness parameter
% sig2: marginal variance
% sig2eps: variance of measurement error

n = length(z);
dist = pdist2(locs, locs);

corr_mat = matern(dist, nu, rho);
L = chol(sig2*corr_mat + sig2eps*speye(n), 'lower');

dist0 = pdist2(pred_locs, locs);
cross_mat = sig2*matern(dist0, nu, rho);

pred_notrend = cross_mat*(L'\(L\z));



if nargout > 1
m = size(pred_locs, 1);
p2 = zeros(m,1);

sig2pred = zeros(m, 1);
for i = 1:m
p2(i) = cross_mat(i,:)*(L'\(L\cross_mat(i,:)'));
sig2pred(i) = sig2 - p2(i);
end

end
