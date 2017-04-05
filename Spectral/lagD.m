function [D,x,w] = lagD(N)
% Laguerre spectral differentiation matrix
[x,w]=gaulag(2,N-1); x=x(:); w=w(:);
c(N-1)=1;
L=[N*(N+1)/2; -x.*LaguerreL(c,3,x).*exp(-x/2)];
X=repmat([0;x], [1,N]);
dX=X-X'+eye(N);
D=(L*(1./L)')./dX-diag(sum(1./dX)+0.5);
w=[2/(N*(N+1)); exp(x-2*log(x)+log(w))];
x=[0; x];
end