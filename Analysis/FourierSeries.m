function [c, L] = FourierSeries(f, a, b, n)
% Computes 2n+1 coefficents of the complex Fourier expansion of f over [a, b].
% Aditionally returns period.
N=2*n;
L=b-a;
x=a+L*(0:N)/N;
h=f(x);
h(1)=(h(1)+h(end))/2;
H=fft(h(1:end-1));
c=1/N*exp(-2i*pi*a/L*(0:N-1)).*H;
c=[c(n+1:end), c(1:n+1)];
end