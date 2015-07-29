function [c, P] = FourierSeries(f, a, b, n)
% Computes 2n+1 coefficents of the complex Fourier expansion of f over [a, b].
% Aditionally returns period.
mid=n+1; n=2*n;
P=b-a;
x=a+P*(0:n)/n;
h=f(x);
h(1)=(h(1)+h(end))/2;
H=fft(h(1:end-1));
c=exp(-2i*pi*a/P*(0:n-1))/n.*H;
c=[c(mid:end), c(1:mid)];
end