function [x, w] = GaussJacobi(a,b,n)
% Returns abscissas and weights for the Gauss-Jacobi n-point quadrature 
% over the interval [-1, 1] using the Golub-Welsch Algorithm.
k=1:n;
c=2*k+a+b;
D=(b^2-a^2)./(c.*(c-2)+(b^2==a^2));
E=2./c.*sqrt(k.*(k+a).*(k+b).*(k+a+b)./(c.^2-1));
[x,V]=trideigs(D, E); 
w0=hypergeom([1,-a],2+b,-1)/(1+b)+hypergeom([1,-b],2+a,-1)/(1+a);
w=w0*V(1,:).^2;
x=x';
end