function [x, w] = GaussHermite(n, mu, sigma)
% Returns abscissas and weights for the Gauss-Hermite n-point quadrature 
% over the interval [-inf, inf] using the Golub-Welsch Algorithm.
E=sqrt((1:n-1)/2);
[x,V]=trideigs(zeros(1,n), E); 
x=(sqrt(2)*sigma)*x'+mu;
w=sqrt(2*pi)*sigma*V(1,:).^2;
end