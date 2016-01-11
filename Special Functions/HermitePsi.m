function [y]=HermitePsi(a, x)
% Evaluates the Hermite function series given by the coefficients a(n) 
n=length(a);
y=(n>1)*a(n); yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=a(k+1)+sqrt(2/(k+1))*x.*y-sqrt((k+1)/(k+2))*yy;
    yy=temp;
end
h0=pi^(-1/4)*exp(-x.^2/2);
y=h0.*(a(1)+sqrt(2)*x.*y-sqrt(1/2)*yy);
end