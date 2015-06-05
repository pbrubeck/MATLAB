function [c, k] = FourierSeries(f, a, b, n)
% Computes n coefficents of the complex Fourier expansion of f over [a, b].
% Adittionally returns frecuency.
P=b-a;
[x, w]=GaussLegendre(a,b,20);
k=(-n:n)*2*pi/P;
c=(f(x).*w/P)*exp(-1i*x(:)*k);
end