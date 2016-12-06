function [D,x,w] = lagD(N)
% Laguerre spectral differentiation matrix
[x,w]=gaulag(1,N-1); x=x(:); w=w(:);
X=repmat([0;x], [1, N]);
c(N-1)=1;
L=[N; -x.*LaguerreL(c,2,x).*exp(-x/2)];
D=(L*(1./L)'-eye(N))./(X-X'+eye(N));
D(1,1)=-N/2;
w=[1/N; exp(x-log(x)+log(w))];
x=[0; x];
end