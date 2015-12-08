function [y]=LaguerreL(c, a, x)
% Evaluates the Laguerre series given by the coefficients c(n)
n=length(c);
b=c(n); bb=zeros(size(x));
for k=n-2:-1:1
    temp=b;
    b=c(k+1)+(2*k+1+a-x)/(k+1).*b-(k+1+a)/(k+2)*bb;
    bb=temp;
end
y=c(1)+(1+a-x).*b-(1+a)/2*bb;
end