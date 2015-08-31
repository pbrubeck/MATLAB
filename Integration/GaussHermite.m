function [x, w] = GaussHermite(n)
% Returns abscissas and weights for the Gauss-Hermite n-point quadrature 
% over the interval [-inf, inf] using the Golub-Welsch Algorithm.
beta=sqrt((1:n-1)/2);
J=full(gallery('tridiag', beta, zeros(1,n), beta));
[V,x]=eig(J,'nobalance','vector');
x=x'; w=sqrt(pi)*V(1,:).^2;
end