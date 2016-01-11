function [x, w] = GaussChebyshev(a, b, n)
% Returns abscissas and weights for the Gauss-Chebyshev n-point quadrature
% over the interval [a, b].
th=pi*(1:2:2*n-1)/(2*n);
x=(b-a)/2*cos(th)+(a+b)/2;
w(1:n)=(b-a)*pi/(2*n);
end