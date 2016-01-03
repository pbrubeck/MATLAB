function [y]=HermiteH(a, x)
% Evaluates the Hermite series given by the coefficients a(n) 
n=length(a);
y=(n>1)*a(n); yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=a(k+1)+2*x.*y-2*(k+1)*yy;
    yy=temp;
end
y=a(1)+2*x.*y-2*yy;
end