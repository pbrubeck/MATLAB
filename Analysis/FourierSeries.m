function [c, k] = FourierSeries(f, a, b, n)
% Computes n coefficents of the complex Fourier expansion of f over [a, b].
% Adittionally returns frecuency.
P=b-a;
h=P/(2*n);
k=2*pi/P*(-n:n);

x=a+h*(0:2*n-1);
w=h*ones(1,2*n);
w(1)=h/2;

s=w.*f(x);
s=fft(s)+h/2*f(b);
J=[s(n+1:end), s(1:n+1)];
c=exp(-1i*a*k)/P.*J;

%[x, w]=GaussChebyshev(a,b,32);
%c=(f(x).*w/P)*exp(-1i*x(:)*k);
end
