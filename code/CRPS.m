	
% Author: Pulong Ma University of Cincinnati
% Date: April 03, 2016 
% Last Modified by: Pulong Ma
% Last Modified Date: April 03, 2016
% Last Modified time: 23:03:08


% Purpose: The continuous ranked probability score for normal distribution



function [crps] = CRPS(x, mu, sig)

%% Input Arguments:
% x: a vector of predicted values
% mu: a vector of predictive mean 
% sig: a vector of predictive standard deviation

xo = (x-mu)./sig;

crps = sig.*(xo.*(2*normcdf(xo)-1) + 2*normpdf(xo) - 1/sqrt(pi)); 	



end
