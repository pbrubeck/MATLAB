function [s] = piecewiseLagrange(x,jumps)
% Correction term due to jumps
m=0:length(jumps)-1;
g=@(xi) Horner((jumps(:)'./factorial(m)), x(:)'-xi)';
s=@(z,xi)  ((z>=xi).*(xi>x)-(z<xi).*(xi<x)).*g(xi);
end