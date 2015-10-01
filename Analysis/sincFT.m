function [u] = sincFT(f, a, b, N)
% Computes the semidiscrete Fourier transform
x0=(b+a)/2;
h=(b-a)/(2*N);
k=(-N:N);
x=k*h-x0;
y=f(x);
u=h/sqrt(2*pi)*((-1).^k(1:end-1)).*(fft(y(1:end-1))+y(end));
u=[u(N+1:-1:1) u(2*N:-1:N+1)];
omega=2*pi/(b-a)*k;

u=u./(-omega.^2); u(N+1)=0;
v=h/sqrt(2*pi)*((-1).^k(1:end-1)).*(fft(u(1:end-1))+u(end));
v=[v(N+1:2*N) v(1:N+1)];

figure(1);
plot(x, [real(v); imag(v)]);
end