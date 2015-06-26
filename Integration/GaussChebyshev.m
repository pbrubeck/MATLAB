function [x, w] = GaussChebyshev(a, b, n)
% Returns abscissas and weights for the Gauss-Chebyshev n-point quadrature
% over the interval [a, b]. Weigths include reciprocal weight function.
th=((1:n)-1/2)*pi/n;
x=((b-a)*cos(th)+a+b)/2;
w=(b-a)*pi/(2*n)*sin(th);
end