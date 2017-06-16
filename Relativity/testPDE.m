function [] = testPDE(m, n)
if nargin<2
    n=m;
end

% Differential operators
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);
A1=Dx*Dx;
A2=Dy*Dy;
A3=4*exp(-xx.^2-yy.^2).*(1-xx.^2-yy.^2);
F=-4*exp(-xx.^2-yy.^2).*(1-exp(-xx.^2-yy.^2)).*(1-xx.^2-yy.^2);

% Boundary condtions
a=[1,1;1,1];
b=[0,0;0,0];
E1=eye(m);
E2=eye(n);
B1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
B2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);

% Right hand sides
f=@(z) exp(-abs(z).^2);
fx=@(z) -2*real(z).*f(z);
fy=@(z) -2*imag(z).*f(z);

w1=[1+1i*y; -1+1i*y];
w2=[x+1i, x-1i];
b1=diag(a(1,:))*f(w1)+diag(b(1,:))*fx(w1);
b2=f(w2)*diag(a(2,:))+fy(w2)*diag(b(2,:));

% Differential operator
opA=@(uu) A1*uu+uu*A2'+A3.*uu;
cA=@(i,j) A1(:,i)+A3(i,j)*(1:m==i)';
rA=@(i,j) A2(:,j)+A3(i,j)*(1:n==j)';

% Preconditioner (Low frequency component)
[P1,P2]=blockprecond(cA,rA,m,n);

% Solver
[green,ps,kd,gb]=elliptic(P1,P2,B1,B2,[1,m],[1,n]);

afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(green(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+green(rhs));

[uu,res,its]=precond(afun,pfun,rhs,uu,20,1e-15);
uu=gb(uu)+ub;

display(its);
display(res);

xq=linspace(-1,1,m);
yq=linspace(-1,1,n);
[xxx,yyy]=ndgrid(xq, yq);
%uuu=interp2(xx',yy',uu',xxx',yyy','splines')';
%uuu=interpcheb2(uu,xxx,yyy);
uuu=interpcheb(interpcheb(uu,xq,1),yq,2);

err=norm(uu-f(xx+1i*yy),'inf');
errint=norm(uuu-f(xxx+1i*yyy),'inf');

figure(1);
surf(xxx,yyy,uuu);
colormap(jet(256));
shading interp;
camlight;

display(err);
display(errint);
end