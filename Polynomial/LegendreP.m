function [y] = LegendreP(a, x)
% Evaluates the Legendre series given by the coeficients a(n)
n=length(a);
b=(n>1)*a(n); bb=zeros(size(x));
for k=n-2:-1:1
    temp=b;
    b=a(k+1)+(2*k+1)/(k+1)*x.*b-(k+1)/(k+2)*bb;
    bb=temp;
end
y=a(1)+x.*b-bb/2;
end