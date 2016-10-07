function [D,t] = fourD(M)
% First derivative matrix with periodic boundary conditions.
dt=2*pi/M;
t=dt*(1:M);
r=[0, 0.5*(-1).^(1:M-1).*cot(dt*(1:M-1)/2)];
D=toeplitz(r,-r);
end

