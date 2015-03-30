function x=fwdSubst(R, b)
n=size(R,1);
x=zeros(n,1);
for i=1:n
    x(i)=(b(i)-R(i,1:i-1)*x(1:i-1))/R(i,i);
end
end