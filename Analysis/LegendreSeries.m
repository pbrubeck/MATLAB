function [c, Q] = LegendreSeries(f, a, b, n)
% Calculates the Legendre expansion of f over the interval [a,b].
% Additionally returns the equivalent polynomial.
[x, w]=GaussLegendre(-1, 1, n+1);
P=Legendre(n);
y1=Horner(P, x);
y2=f(((b-a)*x+b+a)/2);
J=polyInt(polyMult(P,P), -1, 1);
c=(bsxfun(@times, y1, y2)*w(:))./J;
Q=polyComp(c'*P, [a+b -2]/(a-b));
end