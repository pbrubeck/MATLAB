function [c, Q] = ChebyshevSeries(f, a, b, n)
% Calculates the Chebyshev expansion of f over the interval [a,b].
% Uses the discrete cosine transform.
th=((1:n)-1/2)*pi/n;
x=(b-a)*cos(th)+a;
c=sqrt(2/n)*dct(f(x));
c(1)=c(1)/sqrt(2);
T=Chebyshev(n);
Q=polyComp(c*T, [-a 1]/(b-a));
end