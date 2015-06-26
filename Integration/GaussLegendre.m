function [x, w] = GaussLegendre(a, b, n)
% Returns abscissas and weights for the Gauss-Legendre n-point quadrature 
% over the interval [a, b] using the Golub-Welsch Algorithm.

%% Roots of the Legendre polynomial
%P=Legendre(n+1);
%x=polyRoots(P(n+1,:));
%dP=polyD(P(n+1,:));
%w=(b-a)./((1-x.^2).*(Horner(dP, x)).^2);

%% Eigendecomposition of the recurrence matrix
k=(1:n-1);
beta=k./sqrt(4*k.*k-1);
J=full(gallery('tridiag', beta, zeros(1,n), beta));
[V,D]=eig(J);
x=(b-a)/2*diag(D).'+(a+b)/2;
w=(b-a)*V(1,:).^2;
end