function [u_hat] = sincFT(u, L)
% Computes the semidiscrete Fourier transform on [-L/2, L/2)
N=length(u);
h=L/N;
n=-N/2:N/2-1;

u=fftshift(u);
u_hat=h/sqrt(2*pi)*fftshift(fft(u));
omega=2*pi/L*n;

figure(1);
plot(omega, [real(u_hat); imag(u_hat)]);
end