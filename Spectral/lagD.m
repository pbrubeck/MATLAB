function [D,x,w] = lagD(N)
% Laguerre spectral differentiation matrix
[x,w]=GaussLaguerre(0,N); x=x(:);
c=zeros(1,N); c(N)=1;
X=repmat(x,[1, N]);
L=repmat(LaguerreL(c,1,x).*exp(-x/2),[1, N]);
D=(L./L'-eye(N))./(X-X'+eye(N));
end

