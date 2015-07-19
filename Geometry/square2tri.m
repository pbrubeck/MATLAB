function []=square2tri(n)

x=zeros(1,n);
y=zeros(1,n);
for i=0:n-1
    u=i/(n-1);
    for j=0:n-1
        v=j/(n-1);
        m=max(abs(u),abs(v))/(abs(u)+abs(v));
        x(i*n+j+1)=m*u;
        y(i*n+j+1)=m*v;
    end
end
scatter(x,y);
end