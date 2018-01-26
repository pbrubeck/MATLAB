function [s] = piecewiseLagrange(x,jumps)
%PIECEWISELAGRANGE Summary of this function goes here
m=0:length(jumps)-1;
g=@(xi) Horner((jumps(:)'./factorial(m)), x(:)'-xi)';
s=@(z,xi)  ((z>=xi).*(xi>x)-(z<xi).*(xi<x)).*g(xi);
end