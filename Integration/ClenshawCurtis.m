function [x, w] = ClenshawCurtis(a, b, n)
% Returns abscissas and weights for the corresponding Clenshaw-Curtis
% (n+1)-point quadrature over the interval [a, b].
c=1./(1-(0:2:n).^2);
c=[c, c(ceil(n/2):-1:2)];
f=ifft(c, 'symmetric');
w=(b-a)*([f(1)/2, f(2:n), f(1)/2]);
th=(0:n)*pi/n;
x=((b-a)*cos(th)+a+b)/2;
end