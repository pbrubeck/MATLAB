function [A, D, x] = radial(A, D, x)
n=length(x)/2;
x=x(1:n);
A=A(1:n,1:n)+A(1:n,end:-1:n+1);
D=D(1:n,1:n)+D(1:n,end:-1:n+1);
end

