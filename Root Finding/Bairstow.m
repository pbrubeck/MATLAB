function [x] = Bairstow(a)
a=a(:);
u=a([end-1,end-2])/a(end);
n=size(a,1);
D0=eye(n-2);
D1=diag(ones(n-3,1),1);
D2=diag(ones(n-4,1),2);
b=zeros(n,1);
f=zeros(n,1);

err=1;
tol=1e-13;
while err>tol
    A=D0+u(1)*D1+u(2)*D2;
    b(1:end-2)=A\a(3:end);
    f(1:end-2)=A\b(3:end);
    c=a(2)-u'*b(1:2);
    d=a(1)-u(2)*b(1);
    g=b(2)-u'*f(1:2);
    h=b(1)-u(2)*f(1);
    J=[g*u(1)-h,-g;g*u(2),-h];
    u=u-J\[c;d];
    err=max(abs([c;d]));
end
x=(-u(1)+[1;-1]*sqrt(u(1)^2-4*u(2)))/2;

a=b(1:end-2);
switch(min(length(a),4))
    case 4,
        x=[x;Bairstow(a)];
    case 3,
        u=a([2,1])/a(3);
        x=[x;(-u(1)+[1;-1]*sqrt(u(1)^2-4*u(2)))/2];
    case 2,
        x=[x; -a(1)/a(2)];
end
end