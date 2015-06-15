function y = FourierSynth(c, P, x)
% Evaluates Fourier series for a given set of coeficients and period
% using Horner's method.
n=floor(length(c)/2);
th=2i*pi/P*x;
w=exp(th);
z=exp(-n*th);
y=z.*Horner(c, w);
end