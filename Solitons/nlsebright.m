function [u] = nlsebright(c,x0,th0,x,t)
u=sech(x-c*t-x0).*exp(1i*(c*x-(c^2-1)*t/2+th0));
end