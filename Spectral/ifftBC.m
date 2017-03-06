function [u] = ifftBC(v, a, b)
N=length(v);
A=fftBC(eye(N),a,b);
u=A\v;
end