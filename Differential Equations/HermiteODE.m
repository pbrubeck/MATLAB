function [] = HermiteODE(n, m)
%HERMITEODE Summary of this function goes here
%   Detailed explanation goes here
[D,x]=chebD(n);
L=(D-2*diag(x))*D*diag(exp(x.^2));
[V,lambda]=eigs(L(2:end-1,2:end-1),m,'sm');
u=zeros(n,m);
u(2:end-1,:)=V;
plot(x,u); disp(diag(lambda));
end