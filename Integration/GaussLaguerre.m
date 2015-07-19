function [x, w] = GaussLaguerre(n)
% Returns abscissas and weights for the Gauss-Laguerre n-point quadrature
% over the interval [0, inf] using the Golub-Welsch Algorithm.

%% Roots of the Laguerre polynomial
%L=Laguerre(n+2,a);
%x=polyRoots(L(n+1,:));
%w=x./((n+1)*Horner(L(n+2,:), x)).^2;

%% Eigendecomposition of the recurrence matrix
k=(1:n);   alpha=2*k-1; %2*k-1+a;
k=(1:n-1); beta=k;      %sqrt(k.*(k+a));
J=full(gallery('tridiag', beta, alpha, beta));
[V,D]=eig(J);
x=diag(D).';
w=V(1,:).^2;
end