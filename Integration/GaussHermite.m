function [x, w] = GaussHermite(n)
% Returns abscissas and weights for the Gauss-Hermite n-point quadrature 
% over the interval [-inf, inf] using the Golub-Welsch Algorithm.

%% Roots of the Hermite polynomial
%H=Hermite(n+1);
%x=sort(polyRoots(H(n+1,:)));
%w=(2^(n-1)*factorial(n)*sqrt(pi))./(n*Horner(H(n,:), x)).^2;

%% Eigendecomposition of the recurrence matrix
beta=sqrt((1:n-1)/2);
J=full(gallery('tridiag', beta, zeros(1,n), beta));
[V,D]=eig(J);
x=diag(D).';
w=sqrt(pi)*V(1,:).^2;
end