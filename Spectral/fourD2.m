function [D2, t] = fourD2(M)
% Second derivative matrix with periodic boundary conditions.
dt=2*pi/M;
t=dt*(1:M);
D2=toeplitz([-pi^2/(3*dt^2)-1/6, 0.5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
end