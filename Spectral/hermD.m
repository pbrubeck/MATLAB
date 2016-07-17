function [D, x, w] = hermD(N)
% Hermite spectral differentiation matrix
[x,w]=GaussHermite(N+1,0,1/sqrt(2)); x=x(:);
a(N+1)=1;
X=repmat(x,[1, N+1]);
Psi=repmat(HermitePsi(a,x),[1, N+1]);
D=(Psi./Psi'-eye(N+1))./(X-X'+eye(N+1));
end