function [u] = fftBC(v, a, b)
N=length(v);
n=(0:N-1)';
phi=atan2(a,b*n);
v=exp(1i*phi).*v(:);
u=fft(v, 2*N);
u=real(u(1:N));
end