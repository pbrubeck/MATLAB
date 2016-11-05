function [psi,phi]=kink(c,d,x,t)
g=1/sqrt(1-c^2);
psi=4*atan(exp(g*(x-c*t)+d));
phi=-2*g*(c+1).*sech(g*(x-c*t)+d);
end