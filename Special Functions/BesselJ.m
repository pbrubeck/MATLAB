function [J] = BesselJ(a, x)
t=max(x(:));
n=max(32, ceil(2*pi*t));
c=zeros(1,n);
c(1)=exp(-t)/Gamma(a+1);
for i=1:n-1
    c(i+1)=c(i)*t/(i+a);
end
J=LaguerreL(c, a, x.^2/(4*t)).*(x/2).^a;
end