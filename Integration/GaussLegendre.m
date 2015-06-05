function [x, w] = GaussLegendre(a, b, n)
m=(b-a)/2;
c=(a+b)/2;
P=Legendre(n+1);
x=polyRoots(P(n+1,:));
dP=polyD(P(n+1,:));
w=2*m./((1-x.^2).*(Horner(dP, x)).^2);
x=m*x+c;
end
