function [D, x, w] = hermD(N)
% Hermite spectral differentiation matrix
[x,w]=GaussHermite(N,0,1/sqrt(2)); x=x(:);
a(N)=1;
X=repmat(x,[1, N]);
H=repmat(HermitePsi(a,x),[1, N]);
D=(H./H'-eye(N))./(X-X'+eye(N));
end