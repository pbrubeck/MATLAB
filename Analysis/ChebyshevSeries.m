function [c, Q] = ChebyshevSeries(f, a, b, n)
% Calculates the Chebyshev expansion of f over the interval [a,b].
% Uses the discrete cosine transform.
th=((1:n)-0.5)*pi/n;
x=(b-a)/2*cos(th)+(b+a)/2;
c=sqrt(2/n)*dct(f(x));
c(1)=c(1)/sqrt(2);
T=Chebyshev(n);
Q=polyComp(c*T, [a+b -2]/(a-b));
end