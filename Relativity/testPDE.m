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
b1=f([1+1i*y; -1+1i*y]);
b2=f([x+1i, x-1i]);

% Solver
[green,ps,kd,sc,gb]=elliptic(A1,A2,B1,B2,[1,m],[1,n]);

afun=@(uu) sc(uu)+kd(A3).*uu;
pfun=@(rhs) kd(green(rhs));
ub=ps(b1,b2);
rhs=kd(F-A1*ub-ub*A2'-A3.*ub);
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