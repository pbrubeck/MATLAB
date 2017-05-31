function [] = Robin( n )
rd=[1,n];
kd=2:n-1;
[D,x]=chebD(n);
D2=D*D;

alph=[1;1];
beta=[1;-1];
bc=[1;-1];

E=eye(n);
B=diag(alph)*E(rd,:)+diag(beta)*D(rd,:);
G=-B(:,rd)\B(:,kd);
A=D2(kd,kd)+D2(kd,rd)*G;

f=-5*(x.^2-0.25);
rhs=f(kd)-D2(kd,rd)*(B(:,rd)\bc);

u=zeros(n,1);
u(kd)=A\rhs;
u(rd)=G*u(kd)+B(:,rd)\bc;
plot(x,u);
end