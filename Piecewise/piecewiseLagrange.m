function [s1,s2] = piecewiseLagrange(x0,xi,jumps)
% Correction term due to jumps
m=0:length(jumps)-1;
g=Horner(jumps(:)'./factorial(m), x0(:)'-xi)';
s1= -(xi<=x0).*g;
s2=  (xi>=x0).*g;
end