function [x, w] = GaussChebyshev(a, b, n)
% Returns abscissas and weights for the corresponding Gauss-Chebyshev
% n-point quadrature over the interval [a, b].
th=((0:n-1)+1/2)*pi/n;
x=((b-a)*cos(th)+a+b)/2;
w=(b-a)*pi/(2*n)*sin(th);
end