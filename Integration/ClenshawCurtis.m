function [x, w] = ClenshawCurtis(a, b, n)
% Returns abscissas and weights for the corresponding Clenshaw-Curtis
% (n+1)-point quadrature over the interval [a, b].
th=(0:n)*pi/n;
x=(b-a)/2*cos(th)+(b+a)/2;
d=1./(1-(0:2:n).^2);
c=ifft([d, d(ceil(n/2):-1:2)], 'symmetric');
w=(b-a)*[c(1)/2, c(2:n), c(1)/2];
end