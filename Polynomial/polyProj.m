function [c, Q] = polyProj(f, P, a, b)
% Computes the projection of f over P on the interval [a,b]
[x, w]=GaussLegendre(-1,1,20);
y1=Horner(P, x);
y2=arrayfun(f, ((b-a)*x+b+a)/2);
J=polyInt(polyMult(P,P), -1, 1);
c=(bsxfun(@times, y1, y2)*w')./J;
Q=polyComp(c'*P, [a+b -2]/(a-b));
end