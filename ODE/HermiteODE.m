function [lambda] = HermiteODE(n, m)
% Solves ( d^2/dx^2 + x^2 ) Psi = lambda Psi using Hermite spectral methods
[D, x, w]=hermD(n+1);
H=-D*D+diag(x.^2);
[Psi, lambda]=eig(H);
lambda=diag(lambda);
[lambda,order]=sort(lambda);
Psi=Psi(:,order);

% Normalization and interpolation to a finer grid
Psi=normc(Psi, w.*exp(x'.^2));
xx=linspace(x(1), x(end), 4*n); xx=xx(:);
uu=interp1(x,Psi,xx,'spline');
plot(xx,uu(:,m));
title(sprintf('\\lambda_{%d} = %f', m, lambda(m)));
end