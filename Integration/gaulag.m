function [x, w] = gaulag(alpha, n)
% Returns abscissas and weights for the Gauss-Laguerre n-point quadrature
% over the interval [0, inf) using the Golub-Welsch Algorithm.
k=(1:n);
D=2*k-1+alpha;
E=sqrt(k.*(k+alpha));
[x,V]=trideigs(D, E); 
x=x'; w=gamma(alpha+1)*V(1,:).^2;
end