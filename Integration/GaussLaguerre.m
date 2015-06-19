function [x, w] = GaussLaguerre(n)
% Returns abscissas and weights for the corresponding Gauss-Laguerre
% n-point quadrature over the interval [0, inf].
L=Laguerre(n+2);
x=polyRoots(L(n+1,:));
w=x./((n+1)*Horner(L(n+2,:), x)).^2;
end