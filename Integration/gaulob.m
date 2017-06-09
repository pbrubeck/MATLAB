function [x, w] = gaulob(a, b, n)
% Returns abscissas and weights for the Gauss-Lobatto n-point quadrature 
% over the interval [a, b] using the Golub-Welsch Algorithm.
x=zeros(1,n); x([1,n])=[-1,1];
w=zeros(1,n); w([1,n])=(b-a)/(n*(n-1));

k=1:n-3;
E=sqrt((k.*(k+2))./((2*k+1).*(2*k+3)));
[x(2:n-1),V]=trideigs(zeros(1,n-2), E);
w(2:n-1)=2/3*(b-a)*V(1,:).^2./(1-x(2:n-1).^2);
x=(b-a)/2*x+(a+b)/2;
end