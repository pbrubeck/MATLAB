function [u] = sincFT2(f, a, b, N)
% Computes the semidiscrete Fourier transform
x0=(b+a)/2;
h=(b-a)/(2*N);
k=(-N:N);
x=k*h-x0;
[xx, yy]=meshgrid(x);
z=f(xx, yy);
sgn=(-1).^k(1:end-1); 
A=h^2/(2*pi)*(sgn'*sgn);
ii=N+1:-1:1; jj=2*N:-1:N+1;
omega=2*pi/(b-a)*k;
[kx, ky]=meshgrid(omega);

u=A.*fft2(z(1:end-1, 1:end-1));
u=[u(ii, ii), u(ii, jj); u(jj, ii), u(jj, jj)];
% u=-u./(kx.^2+ky.^2); u(N+1,N+1)=0;
% 
% v=(2*pi/(b-a))^2*(sgn'*sgn).*ifft2(u(1:end-1, 1:end-1));
% v=[v(ii, ii), v(ii, jj); v(jj, ii), v(jj, jj)];

colormap(gray(256));
figure(1);
imagesc(omega, omega, real(log(u)));
end