function [] = HermiteODE(n, m)
% Solves ( d^2/dx^2 + x^2 ) Psi = lambda Psi using Hermite spectral methods
[D, x, w]=hermD(n);
[H, L]=eig(-D*D+diag(x.^2), 'vector');
[L,id]=sort(L);
H=H(:,id);

% Normalization and interpolation to a finer grid
H=normc(H, w);
xx=linspace(x(1), x(end), 4*n); xx=xx(:);
uu=interp1(x,H,xx,'spline');
plot(xx,uu(:,m));
end