function [x, w] = GaussChebyshev(a, b, n)
% Returns abscissas and weights for the corresponding Gauss-Chebyshev
% n-point quadrature over the interval [a, b].
m=(b-a)/2;
th=(1:2:2*n-1)*pi/(2*n);
x=m*cos(th)+(a+b)/2;
w=m*pi/n*ones(size(x));
end