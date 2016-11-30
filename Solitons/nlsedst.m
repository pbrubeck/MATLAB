function [x,w,D,Q] = nlsedst(N, dt)
% NonLinear Schrodinger Equation, dst linear propagator
% 1i u_t + 1/2 u_xx = 0

h=sqrt(2*pi/N);
x=h*(-N/2:N/2-1)';
w=h*ones(size(x));
d=1i*sqrt(2*pi/N)*([1:N]');
D=@(u) idst(bsxfun(@times, d, dst(u)));
T=1/2*(2*pi/N)*([1:N]').^2;
q=exp(-1i*dt*T);
Q=@(u) idst(bsxfun(@times, q, dst(u)));
end