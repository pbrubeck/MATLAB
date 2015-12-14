function [y]=ChebT(a, x)
% Evaluates the Chebyshev series given by the coefficients a(n)
n=length(a);
y=(n>1)*a(n); yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=a(k+1)+2*x.*y-yy;
    yy=temp;
end
y=a(1)+x.*y-yy;
end
