function J = ClenshawCurtis(f, a, b, n)
% Returns abscissas and weights for the corresponding Clenshaw-Curtis
% n-point quadrature over the interval [a, b].
th=((0:n-1)+1/2)*pi/n;
x=((b-a)*cos(th)+a+b)/2;
c=dct(f(x));
J=c(1)*sqrt(2);
for k=1:n/2-2
    J=J+2*c(2*k+1)/(1-4*k*k);
end
J=J*(b-a)/sqrt(2*n);
end