function y = FourierSynth(c, P, x)
% Evaluates the Fourier series for a given set of coeficients with period P
% using Horner's method.
n=(length(c)-1)/2;
th=2i*pi/P*x;
w=exp(th);
z=exp(-n*th);
y=z.*Horner(c, w);
end