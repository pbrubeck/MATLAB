function [D,x,w] = lagD(N)
% Laguerre spectral differentiation matrix
[x,w]=GaussLaguerre(1,N-1); x=N/2*x(:);
X=repmat([0;x],[1, N]);
c(N-1)=1;
L=repmat([N; -2/N*x.*LaguerreL(c,2,2/N*x).*exp(-x/N)], [1, N]);
D=N*(L./L'-eye(N))./(X-X'+eye(N));
D(1,1)=-N;
x=[0;x/N];
end