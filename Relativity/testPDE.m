function [] = testPDE(m, n)
% Solves a second-order, self-adjoint differential eqation
if nargin<2
    n=m;
end

% Exact solution
f=@(z) exp(-abs(z*(2+1i)-0.5+0.3i).^2);
fx=@(z) 0*real(z);
fy=@(z) 0*real(z);

% Differential operators
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);

% Eqution coefficients
A11=1+10*x.^2.*y.^2;
A12=x.^2.*y.^2;
A21=x.^2.*y.^2;
A22=1+10*x.^2.*y.^2;
A3=pi^2*sin(x+y);

% Differential operator
opA=@(uu) Dx*(A11.*(Dx*uu)+A12.*(uu*Dy'))+(A21.*(Dx*uu)+A22.*(uu*Dy'))*Dy'+A3.*uu;

% Preconditioner
P1=Dx*diag(mean(A11,2))*Dx+diag(mean(A3,2)-mean(A3(:))/2);
P2=Dy*diag(mean(A22,1))*Dy+diag(mean(A3,1)-mean(A3(:))/2);

% Right-hand side
F=opA(f(xx+1i*yy));

% Boundary condtions
a=[1,1;1,1];
b=[0,0;0,0];
E1=eye(m);
E2=eye(n);
B1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
B2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);

% Right hand sides
yb=[1+1i*y; -1+1i*y];
xb=[x+1i, x-1i];
b1=diag(a(1,:))*f(yb)+diag(b(1,:))*fx(yb);
b2=f(xb)*diag(a(2,:))+fy(xb)*diag(b(2,:));

% Solver
[gf,ps,kd,gb,dL]=elliptic(P1,P2,B1,B2,[1,m],[1,n]);
dL(A11-repmat(mean(A11,2),[1,n]),Dx,E2,Dx,E2);
dL(A22-repmat(mean(A22,1),[m,1]),E1,Dy,E1,Dy);
dL(A12,E1,Dy,Dx,E2);
dL(A21,Dx,E2,E1,Dy);
dL(A3-repmat(mean(A3,2),[1,n])-repmat(mean(A3,1),[m,1])+mean(A3(:)),E1,E2,E1,E2);

afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));

[uu,res,its]=precond(afun,pfun,rhs,uu,100,0e-15);
uu=gb(uu)+ub;

display(its);
display(res);

xq=linspace(-1,1,m);
yq=linspace(-1,1,n);
[xxx,yyy]=ndgrid(xq, yq);
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