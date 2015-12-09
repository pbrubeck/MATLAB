function [y]=HermiteH(a, x)
% Evaluates 
n=length(a);
b=(n>2)*a(n); bb=zeros(size(x));
for k=n-2:-1:1
    temp=b;
    b=a(k+1)+2*x.*b-2*(k+1)*bb;
    bb=temp;
end
y=a(1)+2*x.*b-2*bb;
end