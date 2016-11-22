function [x,w,D,Q] = nlseherm(N, dt)
% NonLinear Schrodinger Equation, Hermite linear propagator
% 1i u_t + 1/2 u_xx = 0

[D,x,w]=hermD(N);
T=-1/2*(D*D);
[V,L]=eig(T(2:end-1,2:end-1),'vector');
L=[0;0;L];
S=zeros(N);
S(2:end-1,3:end)=V;
S(:,1:2)=[0*x+1,x];
Q=S*diag(exp(-1i*dt*L))/S;
end