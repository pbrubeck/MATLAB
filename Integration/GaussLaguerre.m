function [x, w] = GaussLaguerre(n)
% Returns abscissas and weights for the Gauss-Laguerre n-point quadrature
% over the interval [0, inf] using the Golub-Welsch Algorithm.
k=(1:n);   alpha=2*k-1; %2*k-1+a;
k=(1:n-1); beta=k;      %sqrt(k.*(k+a));
J=full(gallery('tridiag', beta, alpha, beta));
[V,D]=eig(J);
x=diag(D).';
w=V(1,:).^2;
end