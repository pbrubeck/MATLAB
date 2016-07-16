function [y] = GegenbauerC(a,alpha,x)
% GegenbauerC(a,alpha,x) Evaluates the Gegenbauer series given by the coeficients a(n)
n=length(a);
y=(n>1)*a(n);
yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=a(k+1)+(2*(k+alpha)/(k+1)*x).*y-(k+2*alpha)/(k+2)*yy;
    yy=temp;
end
y=a(1)+alpha*(2*x.*y-yy);
end