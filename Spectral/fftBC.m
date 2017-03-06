function [u] = fftBC(v, a, b)
N=length(v);
n=(0:N-1)';
phi=atan2(a,b*n);
v1=bsxfun(@times, exp( 1i*phi)/2, v);
v2=bsxfun(@times, exp(-1i*phi)/2, v);
u=fft([v(1,:); v1(2:end,:); zeros(1,size(v,2)); v2(end:-1:2,:)]);
u=u(1:N,:);
end 