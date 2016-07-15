function [y] = JacobiP(c,a,b,x)
% JacobiP(c,a,b,x) Evaluates the Jacobi series given by the coeficients c(n)
n=length(c);
y=(n>1)*c(n);
yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    p=(a+b+2*k+1)/(2*(k+1)*(a+b+k+1))*((a+b+2*k+2)*x+(a^2-b^2)/(a+b+2*k));
    r=-(a+k+1)*(b+k+1)*(a+b+2*k+4)/((k+2)*(a+b+k+2)*(a+b+2*k+2));
    y=c(k+1)+p.*y+r*yy;
    yy=temp;
end
y=c(1)+((a+b+2)*x+a-b).*y/2-(a+1)*(b+1)*(a+b+4)/(2*(a+b+2)^2)*yy;
end