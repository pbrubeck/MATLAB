function [x, w] = GaussLegendre(a, b, n)
% Returns abscissas and weights for the Gauss-Legendre n-point quadrature 
% over the interval [a, b] using the Golub-Welsch Algorithm.
k=(1:n-1);
beta=k./sqrt(4*k.*k-1);
J=full(gallery('tridiag', beta, zeros(1,n), beta));
[V,x]=eig(J,'nobalance','vector');
w=V(1,:).^2;
m=floor(n/2);
x(end:-1:n-m+1)=(x(end:-1:n-m+1)-x(1:n-m))/2;
x(1:n-m)=-x(end:-1:n-m+1);
x=(b-a)/2*x'+(a+b)/2;
w(end:-1:n-m+1)=(w(end:-1:n-m+1)+w(1:n-m))/2;
w(1:n-m)=w(end:-1:n-m+1);
w=(b-a)*w;
end