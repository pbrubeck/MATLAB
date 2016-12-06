function [x, w] = gauleg(a, b, n)
% Returns abscissas and weights for the Gauss-Legendre n-point quadrature 
% over the interval [a, b] using the Golub-Welsch Algorithm.
k=1:n-1;
E=k./sqrt(4*k.*k-1);
[x,V]=trideigs(zeros(1,n), E); 
w=(b-a)*V(1,:).^2;
x=(b-a)/2*x'+(a+b)/2;
end