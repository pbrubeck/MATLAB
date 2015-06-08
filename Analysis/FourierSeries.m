function [c, k] = FourierSeries(f, a, b, n)
% Computes n coefficents of the complex Fourier expansion of f over [a, b].
% Adittionally returns frecuency.
P=b-a;
[x, w]=GaussChebyshev(a,b,40);
k=2*pi/P*(-n:n);
c=(f(x).*w/P)*exp(-1i*x(:)*k);
end