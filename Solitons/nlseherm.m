function [Q,x] = nlseherm(N, dt)
% NonLinear Schrodinger Equation, Hermite linear propagator
% 1i u_t + 1/2 u_xx = 0

[D,x]=hermD(N);
T=-1/2*(D*D);
[V,L]=eig(T(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(2:end-1,3:end)=V;
U(:,1:2)=[0*x+1,x];
Q=U*diag(exp(-1i*dt*L))/U;
end