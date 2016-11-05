function [] = nlsefft(N, init)
% NLSE NonLinear Schrodinger Equation, fft method
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
nframes=1024;
t0=-pi;
tf=pi;

dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

x=sqrt(2*pi/N)*(-N/2:N/2-1)';

% Linear propagator spectrum
% T=-1/2*(D*D); Q=expm(-1i*dt/2*T);
T=1/2*(2*pi/N)*([0:N/2-1, -N/2:-1]').^2;
q=exp(-1i*dt/2*T);

Q=@(u) ifft(q.*fft(u));
nlse(init, Q, x, dt, t0, tf);
end