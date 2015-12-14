function [y]=LaguerreL(c, a, x)
% Evaluates the Laguerre series given by the coefficients c(n)
n=length(c);
y=(n>1)*c(n); yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=c(k+1)+(2*k+1+a-x)/(k+1).*y-(k+1+a)/(k+2)*yy;
    yy=temp;
end
y=c(1)+(1+a-x).*y-(1+a)/2*yy;
end