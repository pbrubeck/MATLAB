function [x, w] = GaussHermite(n)
% Returns abscissas and weights for the Gauss-Hermite n-point quadrature 
% over the interval [-inf, inf] using the Golub-Welsch Algorithm.
beta=sqrt((1:n-1)/2);
[x,V]=trideigs(zeros(1,n), beta); 
x=x'; w=sqrt(pi)*V(1,:).^2;
end