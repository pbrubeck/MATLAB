function [c, k] = FourierSeries(f, a, b, n)
% Computes n coefficents of the complex Fourier expansion of f over [a, b].
% Adittionally returns frecuency.
P=b-a;
k=2*pi/P*(-n:n);
x=a+P*(0:2*n-1)/(2*n);

h=f(x);
h(1)=h(1)/2;
H=fft(h)+f(b)/2;
H=[H(n+1:end), H(1:n+1)];
c=exp(-1i*a*k)/(2*n).*H;

%[x, w]=GaussChebyshev(a,b,32);
%c=(f(x).*w/P)*exp(-1i*x(:)*k);
end
