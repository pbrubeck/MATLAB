function [x, w] = gaujac(a,b,n)
% Returns abscissas and weights for the Gauss-Jacobi n-point quadrature 
% over the interval [a, b] using the Golub-Welsch Algorithm.
k=1:n;
c=2*k+a+b;
D=(b^2-a^2)./(c.*(c-2)+(b^2==a^2));
E=2./c.*sqrt(k.*(k+a).*(k+b).*(k+a+b)./(c.^2-1));
[x,V]=trideigs(D, E); 
w0=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
w=w0*V(1,:).^2;
x=x';
end