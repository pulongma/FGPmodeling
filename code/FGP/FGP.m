% Author: Pulong Ma
% Date: August 31, 2015 
% Last Modified by: Pulong Ma	
% Last Modified Date: February 11, 2016
% Last Modified time: 22:52:48


% Purpose: Spatial prediction in SRE + CAR model for all the 
% locations (including observed and unobserved)

function [pred, sig2SIK] = FGP(data, pred_locs, Xp, K, tau2, ...
		gamma, H, A, basis_center, threshold, sig2eps, Delta)
% observed data
lat = data(:,1);
lon = data(:,2);
z = data(:,3);
Xo = data(:,4:end);



% coordinates of prediction locations
lat_pred = pred_locs(:,1);
lon_pred = pred_locs(:,2); 


n = size(z,1);  % number of observations
m = length(lon_pred); % number of prediction locations
N = size(H, 1);
Veps = sig2eps*speye(n);

if nargin < 12 % homogeneous CAR 
	Delta = speye(N);
	DeltaInv = speye(N);
else % weighted CAR
	DeltaInv = spdiags(diag(Delta).^(-1), 0, N, N);
end	

% BFs evaluated at observed locations
So = bisquare([lat, lon], basis_center, threshold);

% BFs evaluated at prediction locations
Sp = bisquare([lat_pred, lon_pred], basis_center, threshold);
r = size(So, 2);

% link matrix
Ao = A{1};
Ap = A{2};

% help terms
VInv = spdiags((diag(Veps)).^(-1), 0, n, n);
AVA = Ao'*VInv*Ao;
AVz = Ao'*VInv*z;
AVS = Ao'*(VInv*So);

VInvX = VInv*Xo;
AVX = Ao.'*VInvX;


%%%%%%%%%%  smooth the data  %%%%%%%%%%%%%%%
Q = (DeltaInv - gamma*H)/tau2;
mid = Q + AVA;
[L1, ~, s1] = lchol(mid); clear mid;
mid4S(s1, :) = L1'\(L1\AVS(s1, :)); clear AVS;
DSo = VInv*So - VInv*(Ao*mid4S); clear mid4S;
SDS = So'*DSo;
temp = K\speye(r) + SDS;
%mid3z(s1, 1) = L1'\(L1\AVz(s1,1));
mid3z(s1, 1) = cs_ltsolve(L1, L1\AVz(s1,1));
Dz = VInv*z - VInv*(Ao*mid3z);
SigInvZ = Dz - DSo*(temp\(So'*Dz));


% predict smoothed residuals
part1 = Sp*(K*(So'*SigInvZ)); 
[Lq, ~, sq] = lchol(Q);
ASigInvZ = Ao'*SigInvZ;
%QASigInvZ(sq, 1) = Lq'\(Lq\ASigInvZ(sq,1));
QASigInvZ(sq, 1) = cs_ltsolve(Lq, Lq\ASigInvZ(sq,1));

part2 = Ap*QASigInvZ;
pred = part1 + part2;




%%%%%%%%  calculate the SIK variance  %%%%%%%%
if nargout > 1  % upon request
t0 = tic;	

SpK = Sp*K;
% part: uncertainties in estimating regression coefficients
numcoeff = size(Xo,2);
mid3X(s1, :) = L1'\(L1\AVX(s1,:));
DX = VInvX - VInv*Ao*mid3X;
SigInvX = DX - DSo*(temp\(So'*DX));
XSigX_chol = chol(Xo'*SigInvX, 'lower') \ eye(numcoeff);

tempX1 = SpK*(So'*SigInvX);
ASigInvX = Ao'*SigInvX;
QASigX(sq, :) = Lq'\(Lq\ASigInvX(sq,:));

tempX2 = Ap*QASigX;

CpSigInvX = tempX1 + tempX2;

tempX = (Xp - CpSigInvX)*XSigX_chol;

sig2_beta = sum(tempX.*tempX, 2);





% part1: Sp*K*So'*D*So*K*Sp - SpK*So'*D*So*temp^(-1)*SDS*SpK'
p0 = sum(SpK.*Sp, 2); Sp=[];
tempSDS = temp\SDS;
p1 = sum((SpK*SDS*(speye(r)-tempSDS)).*SpK, 2);

% part2: SpK*So'*D*Ao*Q^(-1)*Ap' + Ap*Q^(-1)*Ao'*DSo*SpK'
%  -[SpK*SDS*tempSo'*D*Ao*Q^(-1)*Ap'+ Ap*Q^(-1)*Ao'*DSo*tempSDS*SpK']
ADS = Ao'*DSo; clear DSo;
QInvADS(sq, :) = Lq'\(Lq\ADS(sq, :)); clear ADS;
ApQInvADS = Ap*QInvADS; clear QInvADS;
p2 = 2*sum((ApQInvADS*(speye(r)-tempSDS)).*SpK, 2);
clear SpK;

% part4: Ap*Q^(-1)*Ao'*D*So*temp^(-1)*So'*D*Ao*Q^(-1)*Ap';
p4 = sum((ApQInvADS/temp).*ApQInvADS, 2);
clear ApQInvADS;

% part3: Ap*Q^(-1)*Ao'*D*Ao*Q^(-1)*Ap'
p3 = zeros(m,1);
sig2xi = zeros(m, 1);
%for i = 1:m
%	QAp(sq, 1) = cs_ltsolve(Lq, Lq\full(Ap(i,sq)'));
%	sig2xi(i) = Ap(i,:)*QAp;
%	AVAQAp = AVA*QAp;
%	mid4QAp(s1, 1) = cs_ltsolve(L1, L1\AVAQAp(s1,1));
%	DAQAp = VInv*(Ao*(QAp-mid4QAp));
%	p3(i) = (Ao*QAp)'*DAQAp;
%end

if m > 1
	if m*N*8/10^9<5
		BLOCKSIZE = 200;
	else
		BLOCKSIZE = 100;
	end
	counter = 0;
	while(counter < m)
		ind = [(counter+1):min(counter+BLOCKSIZE,m)].';
		QAp(sq, :) = Lq.'\(Lq\Ap(ind,sq).');
		sig2xi(ind) = sum(Ap(ind,:).*QAp.', 2);
		AVAQAp = AVA*QAp;
		mid4QAp(s1, :) = L1.'\(L1\AVAQAp(s1, :));
		DAQAp = VInv*(Ao*(QAp-mid4QAp));
		p3(ind) = sum((Ao*QAp).*DAQAp, 1).';
		counter = counter + BLOCKSIZE;
		clear QAp mid4QAp;
	end
else % if there is 
	QAp(sq, 1) = cs_ltsolve(Lq, Lq\full(Ap(sq)'));
	sig2xi = Ap*QAp;
	AVAQAp = AVA*QAp;
	mid4QAp(s1, 1) = cs_ltsolve(L1, L1\AVAQAp(s1));
	DAQAp = VInv*(Ao*(QAp-mid4QAp));
	p3 = (Ao*QAp)'*DAQAp;
end

% putting all pieces together
sig2SIK = p0+sig2xi - (p1+p2+p3-p4) + sig2_beta;
t=toc(t0);
% disp(strcat('time to compute kriging variance:', num2str(t/60), 'mins'))
end



end % end main function


