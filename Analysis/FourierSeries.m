function [c, P] = FourierSeries(f, a, b, n)
% Computes 2n+1 coefficents of the complex Fourier expansion of f over [a, b].
% Adittionally returns period.
n=2*n;
P=b-a;
x=a+P*(0:n)/n;
h=f(x);
fa=h(1);
fb=h(end);
h(1)=(fa+fb)/2;
H=fft(h(1:end-1));
c=exp(-2i*pi*a/P*(0:n-1))/n.*H;
c=[c(n/2+1:end), c(1:n/2+1)];
end
