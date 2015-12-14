function [y]=Clenshaw(A, B, C, a, x)
% Evaluates the recursion P_n(x)=(A_n*x+B_n)*P_{n-1}(x)+C_n*P_{n-2}(x)
% Assumes P_0(x)=1 and P_1(x)=x
n=length(a);
A(1:n)=A; B(1:n)=B; C(1:n)=C;
y=(n>2)*a(n); yy=zeros(size(x));
for k=n-2:-1:1
    temp=y;
    y=a(k+1)+(A(k)*x+B(k)).*y+C(k+1)*yy;
    yy=temp;
end
y=a(1)+x.*y+C(1)*yy;
end