function [u,v]=kink(c,d,x,t)
g=1/sqrt(1-c^2);
u=4*atan(exp(g*(x-c*t)+d));
v=-2*g*(c+1).*sech(g*(x-c*t)+d);
end