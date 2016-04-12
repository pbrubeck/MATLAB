function [a,p] = Frobenius(b,c)
% Solves the general 2nd order linear ODE through Frobenius Method
% y'' + g1(x)/(x-x0)*y' + g2(x)/(x-x0)^2*y=0
% Where g1(x)=sum(b_n (x-x0)^n) and g2(x)=sum(c_n (x-x0)^n)
% Outputs the solution with the form y(x)=sum(a_n (x-x0)^(n+p))
% Important: note that a(1)=a_0, the same holds for vectors b and c
b=b(:); 
c=c(:);
N=length(b);
a=zeros(1,N);
a(1)=1;
a(2)=0;
p=(b(1)-1+sqrt((1-b(1))^2-4*c(1)))/2;
for n=1:N-1
    a(n+1)=-a(n:-1:1)*((n+p-(1:n)').*b(2:n+1)+c(2:n+1))/((n+p)*(n+p-1+b(1))+c(1));
end
end