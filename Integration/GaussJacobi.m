function [x, w] = GaussJacobi(alpha,beta,n)
% Returns abscissas and weights for the Gauss-Jacobi n-point quadrature 
% over the interval [a, b] using the Golub-Welsch Algorithm.
a=-1;
b=1;
k=(1:n);
c=2*k+alpha+beta;
D=(beta^2-alpha^2)./(c.*(c-2)+(beta^2==alpha^2));
E=2./c.*sqrt(k.*(k+alpha).*(k+beta).*(k+alpha+beta)./(c.^2-1));
[x,V]=trideigs(D, E); 
w=(b-a)*V(1,:).^2;
x=(b-a)/2*x'+(a+b)/2;
end