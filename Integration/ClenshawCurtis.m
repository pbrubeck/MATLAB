function [x, w] = ClenshawCurtis(a, b, n)
% Returns abscissas and weights for the corresponding Clenshaw-Curtis
% n-point quadrature over the interval [a, b].
th=(0:n-1)*pi/(n-1);
x=(b-a)/2*cos(th)+(b+a)/2;
d=1./(1-(0:2:n-1).^2);
c=ifft([d, d(ceil((n-1)/2):-1:2)], 'symmetric');
w=(b-a)*[c(1)/2, c(2:n-1), c(1)/2];
end