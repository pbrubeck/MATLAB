function [u]=nlsedark(k,phi,x0,th0,x,t)
a=cos(phi);
v=sin(phi)+k;
u=(a*tanh(a*(x-v*t-x0))+1i*(v-k)).*exp(1i*(k*x-(k^2/2+a^2+(v-k)^2)*t+th0));
end