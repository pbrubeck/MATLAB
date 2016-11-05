function [Q,x] = nlsefft(N, dt)
% NonLinear Schrodinger Equation, fft linear propagator
% 1i u_t + 1/2 u_xx = 0

x=sqrt(2*pi/N)*(-N/2:N/2-1)';
T=1/2*(2*pi/N)*([0:N/2-1, -N/2:-1]').^2;
q=exp(-1i*dt*T);
Q=@(u) ifft(q.*fft(u));
end