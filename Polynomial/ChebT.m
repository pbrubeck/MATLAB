function [y]=ChebT(a, x)
% Evaluates the Chebyshev series given by the coefficients a(n)
n=length(a);
b=(n>1)*a(n); bb=zeros(size(x));
for k=n-2:-1:1
    temp=b;
    b=a(k+1)+2*x.*b-bb;
    bb=temp;
end
y=a(1)+x.*b-bb;
end
