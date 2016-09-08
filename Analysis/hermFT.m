function [F] = hermFT(x,w)
% Computes Fourier Transform operator given the Gauss-Hermite nodes and weights
F=1/sqrt(2*pi)*exp(-1i*(x(:)*x(:)'))*diag(w);
end