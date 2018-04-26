function [s1,s2] = piecewiseLagrange(x0,xi,jumps)
% Correction Lagrange coefficients due to jumps
m=0:size(jumps,1)-1;
g=Horner((jumps.')./factorial(m), x0(:).'-xi).';
s1 = -(xi<=x0(:)).*g;
s2 =  (xi>=x0(:)).*g;
end