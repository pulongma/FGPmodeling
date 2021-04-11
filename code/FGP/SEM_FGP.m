	
% Author: Pulong Ma
% Date: June 02, 2016 
% Last Modified by: Pulong Ma
% Last Modified Date: June 08, 2016
% Last Modified time: 15:19:29


% Purpose: EM algorithm is used to estimate parameters in 
% SRE + CAR model. eta and xi are considered to be missing data. 


function [K_em,tau2_em,gamma_em,beta_em,sig2eps_em,T,epsilon,time] = SEM_SIK(S,...
	     A,Z,X,Hall,sig2eps,Delta,maxit,avgtol,opts)

%% Input Arguments:
% S: is an nxr basis function matrix.
% Z: is an nx1 data vector.
% X: is an nxp matrix of covariates.
% Hall: is a struct array containing the NxN proximity matrix 
% on discretized domain and smallest and largest eigenvalues 
% sig2eps: is the measurement error. When it is not given, then the function 
% can estiamte this measurement error; otherwise this algorithm will use 
% the true value of sig2eps given by the user.
% Delta: is the NxN diagonal matrix
% maxit: is the maximum number of iterations in MCEM algorithm.
% avgtol: the tolerance error when estimating parameters.
% opts: is an optimzation object specifying the optimzation algorithm.


H = Hall.H; % NxN proximity matrix
n = length(Z); % number of observed data
r = size(S,2); % number of basis function
N = size(H,1); % number of area units in discretized domain
p = size(X,2);

% default values for the last two parameters
if nargin<9, avgtol=1e-4; end
if nargin<8, maxit=250; end
if nargin<7 % homogeneous CAR 
	Delta = speye(N);
	DeltaInv = speye(N);
	large = Hall.large;
	small = Hall.small;                                                
	lbg = 1/small + 1e-3*abs(1/small);
	ubg = 1/large - 1e-3*abs(1/large);
	gamma = (lbg + ubg)/2;
else % weighted CAR
	DeltaInv = spdiags(diag(Delta).^(-1), 0, N, N);
	lbg = -1+1e-3;
	ubg = 1-1e-3;
	gamma = 0.16;
end	

% initial values
varest = var(Z, 1);
K_old = 0.95*varest*eye(r);
tau2 = 0.02*varest;
t = 1;
done = 0;
neg2loglikelihood = NaN;
beta_old = regress(Z,X);



update = 0;
lb = lbg;
ub = ubg;

if nargin < 10
	opts = optimoptions(@fmincon,'TolX', 1e-3, 'Algorithm', ...
		'active-set','Display','off');
end


fun = @neg2Q_gamma;



L = 1; % SEM
beta_new = zeros(p);
K_new = zeros(r,r);


%% if measurement-error variance is provided
if ~isempty(sig2eps)  
% help term
Veps = sig2eps*speye(n);
VInv = spdiags((diag(Veps)).^(-1), 0, n, n);
AVA = A'*VInv*A;
AVS = A'*(VInv*S);
XVXInv = (X'*VInv*X)\eye(p);
XVInv = X'*VInv;

ts=tic;

while done == 0


%% compute posterior mean and variance for eta
ztilde = Z - X*beta_old;
AVztilde = A'*(VInv*ztilde);
Q = (DeltaInv-gamma(t)*H)/tau2(t);
mid = Q + AVA;
[Lc, ~, sc] = lchol(mid);
%[Lc, ~, sc] = chol(mid, 'lower', 'vector');
mid4S(sc, :) = Lc'\(Lc\AVS(sc, :));
DS = VInv*S - VInv*(A*mid4S);
SDS = S'*DS; clear mid4S;
temp = K_old\speye(r) + SDS;
mid3ztilde(sc, 1) = cs_ltsolve(Lc, Lc\AVztilde(sc,1)); 
%mid3ztilde(sc, 1) = Lc'\(Lc\AVztilde(sc,1));
Dztilde = VInv*ztilde - VInv*(A*mid3ztilde);

% posterior mean and variance of eta given Z
muEta = K_old*(S'*Dztilde) - K_old*(SDS*(temp\(S'*Dztilde)));
SigEta = K_old - K_old*SDS*K_old' + K_old*SDS*(temp\SDS)*K_old';
% update K 
K_new = SigEta + muEta*muEta';


%% conditional simulation for xi given Z
[Lq,~,sq] = lchol(Q);
%[Lq,~,sq] = chol(Q, 'lower', 'vector');
b = randn(N,1);
xi_NS(sq,1) = cs_ltsolve(Lq,b(sq,1));
%xi_NS(sq,1) = Lq'\b(sq,1);
z = ztilde - S*muEta - A*xi_NS;
Vz = VInv*z;
AVz = A'*Vz;
CInvAvz(sc,1) = cs_ltsolve(Lc, Lc\AVz(sc,1));
%CInvAvz(sc,1) = Lc'\(Lc\AVz(sc,1));
Dz = Vz - VInv*(A*CInvAvz);
SigInvz = Dz - DS*(temp\(DS'*z));
ASigInvz = A'*SigInvz;
QASigInvz(sq,1) = cs_ltsolve(Lq, Lq\ASigInvz(sq,1));
%QASigInvz(sq,1) = Lq'\(Lq\ASigInvz(sq,1));
xi_CS = xi_NS + QASigInvz;

% update beta 
beta_new = XVXInv*(XVInv*(Z-S*muEta-A*xi_CS));


% update tau2
xiHxi = xi_CS'*H*xi_CS;
tau2(t+1) = (xi_CS'*DeltaInv*xi_CS - gamma(t)*xiHxi)/N;

% update gamma

% find numerical solution for gamma
if update == 0
	theta0 = gamma(t);
	theta = fmincon(fun, theta0, [], [], [], [], lb, ub, [], opts);
	gamma(t+1) = theta;
elseif update == 1
	gamma(t+1) = gamma(t);
end

% check convergence 

if abs(gamma(t+1)-gamma(t)) < 1e-3 & t>5
	update = 1;
end

diff = sum(sum((K_new-K_old).^2, 1), 2) + (tau2(t+1)-tau2(t))^2 ...
		+ (gamma(t+1)-gamma(t))^2 + sum((beta_new-beta_old).^2);
if diff < min(avgtol*r^2,1)
	done = 1;
end


if t > maxit
	done = 1;
	% disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
end


beta_old = beta_new;
K_old = K_new;
t = t+1;
%disp(strcat(' t= ', num2str(t)));


end  % end while

te = toc(ts);

sig2eps_em = sig2eps;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if measurement-error variance is not provided

if isempty(sig2eps)  
sig2eps(1) = 0.01*varest;

% help term
Veps_temp = speye(n);
VInv_temp = spdiags((diag(Veps_temp)).^(-1), 0, n, n);
AVA_temp = A'*VInv_temp*A;
AVS_temp = A'*(VInv_temp*S);
XVXInv_temp = (X'*VInv_temp*X)\eye(p);
XVInv_temp = X'*VInv_temp;


ts=tic;

while done == 0

% help term
Veps = sig2eps(t)*Veps_temp;
VInv = (1/sig2eps(t))*VInv_temp;
AVA = (1/sig2eps(t))*AVA_temp;
AVS = (1/sig2eps(t))*AVS_temp;
XVXInv = sig2eps(t)*XVXInv_temp;
XVInv = (1/sig2eps(t))*XVInv_temp;

%% compute posterior mean and variance for eta
ztilde = Z - X*beta_old;
AVztilde = A'*(VInv*ztilde);
Q = (DeltaInv-gamma(t)*H)/tau2(t);
mid = Q + AVA;
[Lc, ~, sc] = lchol(mid);
%[Lc, ~, sc] = chol(mid, 'lower', 'vector');
mid4S(sc, :) = Lc'\(Lc\AVS(sc, :));
DS = VInv*S - VInv*(A*mid4S);
SDS = S'*DS; clear mid4S;
temp = K_old\speye(r) + SDS;
mid3ztilde(sc, 1) = cs_ltsolve(Lc, Lc\AVztilde(sc,1)); 
%mid3ztilde(sc, 1) = Lc'\(Lc\AVztilde(sc,1));
Dztilde = VInv*ztilde - VInv*(A*mid3ztilde);

% posterior mean and variance of eta given Z
muEta = K_old*(S'*Dztilde) - K_old*(SDS*(temp\(S'*Dztilde)));
SigEta = K_old - K_old*SDS*K_old' + K_old*SDS*(temp\SDS)*K_old';
% update K 
K_new = SigEta + muEta*muEta';


%% conditional simulation for xi given Z
[Lq,~,sq] = lchol(Q);
%[Lq,~,sq] = chol(Q, 'lower', 'vector');
b = randn(N,1);
xi_NS(sq,1) = cs_ltsolve(Lq,b(sq,1));
%xi_NS(sq,1) = Lq'\b(sq,1);
z = ztilde - S*muEta - A*xi_NS;
Vz = VInv*z;
AVz = A'*Vz;
CInvAvz(sc,1) = cs_ltsolve(Lc, Lc\AVz(sc,1));
%CInvAvz(sc,1) = Lc'\(Lc\AVz(sc,1));
Dz = Vz - VInv*(A*CInvAvz);
SigInvz = Dz - DS*(temp\(DS'*z));
ASigInvz = A'*SigInvz;
QASigInvz(sq,1) = cs_ltsolve(Lq, Lq\ASigInvz(sq,1));
%QASigInvz(sq,1) = Lq'\(Lq\ASigInvz(sq,1));
xi_CS = xi_NS + QASigInvz;

% update beta 
Ztemp = Z-S*muEta-A*xi_CS;
beta_new = XVXInv*(XVInv*Ztemp);

% update sig2eps
sig2eps(t+1) = (Ztemp - X*beta_new)'*VInv_temp*(Ztemp - X*beta_new)/n;

% update tau2
xiHxi = xi_CS'*H*xi_CS;
tau2(t+1) = (xi_CS'*DeltaInv*xi_CS - gamma(t)*xiHxi)/N;

% update gamma

% find numerical solution for gamma
if update == 0
	theta0 = gamma(t);
	theta = fmincon(fun, theta0, [], [], [], [], lb, ub, [], opts);
	gamma(t+1) = theta;
elseif update == 1
	gamma(t+1) = gamma(t);
end

% check convergence 

if abs(gamma(t+1)-gamma(t)) < 1e-4 & t>5
	update = 1;
end

diff = sum(sum((K_new-K_old).^2, 1), 2) + (tau2(t+1)-tau2(t))^2 ...
		+ (gamma(t+1)-gamma(t))^2 + sum((beta_new-beta_old).^2) ...
		+ (sig2eps(t+1)-sig2eps(t))^2;

if diff < min(avgtol*r^2,1)
	done = 1;
end


if t > maxit
	done = 1;
	% disp(strcat('Algorithm did not converge after ', num2str(maxit),' iterations')); 
end


beta_old = beta_new;
K_old = K_new;

t = t+1;
%disp(strcat(' t= ', num2str(t)));


end  % end while

sig2eps_em = sig2eps(t);

end % end if

te = toc(ts);




% check positive definiteness of K
E = eig(K_new);
if nnz(E<=0) > 0
	disp('Error: K is NOT Positive definite!')
end

K_em = K_new;
tau2_em = tau2(t);
gamma_em = gamma(t);

beta_em = beta_new;


if nargout > 5
	T = t;
end

if nargout > 6
	epsilon = diff;
end	

if nargout > 7
	time = te;
end


function f = neg2Q_gamma(theta)
	Q = (DeltaInv - theta*H);
	[Lq,~,sq] = lchol(Q);
	%[Lq,~,sq] = chol(Q, 'lower', 'vector');
	f = -2*sum(log(diag(Lq))) - theta*xiHxi/tau2(t);
end





end






