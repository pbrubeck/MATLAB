function [u] = sincFT2(f, N)
% Computes the semidiscrete 2D Fourier transform
L=sqrt(2*pi*N);
h=L/N;
n=-N/2:N/2-1;
x=h*n;
k=2*pi/L*n;
[xx, yy]=meshgrid(x);
[kx, ky]=meshgrid(k);

z=f(xx, yy);
A=h^2/(2*pi);
u=fftshift(fft2(fftshift(z))*A); 
%uu=-u./(kx.^2+ky.^2); uu(N/2+1,N/2+1)=0;
v=fftshift(ifft2(fftshift(z))/A);


figure(1); colormap(jet(256));
subplot(1,2,1); imagesc(x,x,real(z)); colorbar();
subplot(1,2,2); imagesc(k,k,real(u)); colorbar();
end