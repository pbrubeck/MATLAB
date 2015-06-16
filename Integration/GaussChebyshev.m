function [x, w] = GaussChebyshev(a, b, n)
% Returns abscissas and weights for the corresponding Gauss-Chebyshev
% n-point quadrature over the interval [a, b].
m=(b-a)/2;
c=(a+b)/2;
th=(1:2:2*n-1)*pi/(2*n);
u=cos(th);
x=m*u+c;
w=m*pi/n*sqrt(1-u.^2);
end
