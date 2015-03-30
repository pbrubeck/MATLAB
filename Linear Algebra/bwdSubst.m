function x=bwdSubst(R, b)
n=size(R,1);
x=zeros(n,1);
for i=n:-1:1
    x(i)=(b(i)-R(i,i+1:end)*x(i+1:end))/R(i,i);
end
end