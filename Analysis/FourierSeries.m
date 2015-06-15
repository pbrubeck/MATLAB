function [c, P] = FourierSeries(f, a, b, n)
% Computes n coefficents of the complex Fourier expansion of f over [a, b].
% Adittionally returns frecuency.
P=b-a;
x=a+P*(0:2*n-1)/(2*n);

h=f(x);
h(1)=(f(a)+f(b))/2;
H=fft(h);
H=[H(n+1:end), H(1:n+1)];
c=exp(-2i*pi*(-n:n)*a/P)/(2*n).*H;
end