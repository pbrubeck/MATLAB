function [x,w,D,Q] = nlsefft(N, dt)
% NonLinear Schrodinger Equation, fft linear propagator
% 1i u_t + 1/2 u_xx = 0

h=sqrt(2*pi/N);
x=h*(-N/2:N/2-1)';
w=h*ones(size(x));
d=1i*sqrt(2*pi/N)*([0:N/2-1, -N/2:-1]');
D=@(u) ifft(bsxfun(@times, d, fft(u)));
T=1/2*(2*pi/N)*([0:N/2-1, -N/2:-1]').^2;
q=exp(-1i*dt*T);
Q=@(u) ifft(bsxfun(@times, q, fft(u)));
end