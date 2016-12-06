function [D, x, w] = hermD(N)
% Hermite spectral differentiation matrix
[x,w]=gauherm(N,0,1/sqrt(2));
x=x(:);
a(N)=1;
h=HermitePsi(a,x);
X=repmat(x,[1, N]);
H=repmat(h,[1, N]);
D=(H./H'-eye(N))./(X-X'+eye(N));
w=1./(N*h.^2);
end