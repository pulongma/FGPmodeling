	
% Author: Pulong Ma
% Date: August 20, 2015 
% Last Modified by: Pulong Ma
% Last Modified Date: August 20, 2015
% Last Modified time: 13:49:24


% Purpose: Calculate Matern correlation function 



function [mat_corr]=matern(dists,smooth,range)

mat_corr=matern_corr(dists/range,smooth);


%%%%%%    the matern correlation function   %%%%%%%
% this parameterization coincides with the double exponential at nu -> Inf
% h is the distance divided by the range parameter
function R=matern_corr(h,nu)

if nu == Inf
	R = exp(-h.^2);
else	
    tau=2*sqrt(nu).*h;
    R = tau.^nu.*besselk(nu,tau)./(2.^(nu-1).*gamma(nu));
    R(h==0)=1;

end

end



end
