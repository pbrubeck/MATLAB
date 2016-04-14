function [df] = fdmconv(f, a, b)
n=length(f);
h=(b-a)/(n-1);
D=[-1,8,0,-8,1]/(12*h);
df=conv(D,f);
df=df(3:end-2);
A=[-3,4,-1;-1,0,1]/(2*h);
B=[-1,0,1;1,-4,3]/(2*h);
df(1:2)=A*f(1:3);
df(n-1:n)=B*f(n-2:n);
end