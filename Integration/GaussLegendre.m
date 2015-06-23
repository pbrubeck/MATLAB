function [x, w] = GaussLegendre(a, b, n)
% Returns abscissas and weights for the corresponding Gauss-Legendre
% n-point quadrature over the interval [a, b].
P=Legendre(n+1);
x=polyRoots(P(n+1,:));
dP=polyD(P(n+1,:));
w=(b-a)./((1-x.^2).*(Horner(dP, x)).^2);
x=((b-a)*x+a+b)/2;
end