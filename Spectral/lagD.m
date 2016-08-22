function [D,x,w] = lagD(N)
% Laguerre spectral differentiation matrix
[x,w]=GaussLaguerre(1,N-1); x=x(:);
X=repmat([0;x], [1, N]);
c(N-1)=1;
L=repmat([N; -x.*LaguerreL(c,2,x).*exp(-x/2)], [1, N]);
D=(L./L'-eye(N))./(X-X'+eye(N));
D(1,1)=-N/2;
x=[0;x];
w=[0,w];
end