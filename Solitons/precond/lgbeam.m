function [uu] = lgbeam(rr,th,p,l,omega)
a=sqrt((omega*factorial(p))/(pi*factorial(p+abs(l))));
c=[zeros(1,p);a];
xx=omega*rr.^2;
uu=(xx.^abs(l/2)).*LaguerreL(c,l,xx).*exp(-xx/2+1i*l*th);
end