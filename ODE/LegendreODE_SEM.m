function [] = LegendreODE_SEM(k,p,m)
mask=(mod(0:(p-1)*(k-1),p-1)==0)';

% Element matrices
[D0,x0]=legD(p);
VL=VandermondeGeg(1/2, x0);
VG=VandermondeGeg(3/2, x0);
M0=inv(VL*VL');
S0=VG\D0;
K0=S0'*S0;

% Element boundaries
xe=linspace(-1,1,k+1)';
J0=diff(xe)/2; % Jacobian
x1=kron(J0, x0)+kron((xe(1:end-1)+xe(2:end))/2, ones(p,1));
x1(p+1:p:end)=[];

% Spectral element assembly
J1=kron(J0,ones(p-1,1));
J1=J1(1:numel(mask));
M1=conv2(diag(mask.*J1), M0);
K1=conv2(diag(mask./J1), K0);

[P,L]=eig(K1,M1, 'vector');
[L,id]=sort(L);
P=P(:,id);

% Normalization
P=bsxfun(@rdivide, P(:,m), P(end,m));
plot(acos(x1),P); xlim([0,pi]);
disp(L(m));
end