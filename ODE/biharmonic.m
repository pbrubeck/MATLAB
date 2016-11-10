function [] = biharmonic( n )
% solves u_xxxx = f(x) subject to u(-1)=u(1)=u'(-1)=u'(1)=0

[D,x]=chebD(n);
L=(diag(1-x.^2)*D^4-8*diag(x)*D^3-12*D^2)*diag(1./(1-x.^2));
f=exp(x);
u=[0; L(2:end-1,2:end-1)\f(2:end-1); 0];

plot(x,u);
end

