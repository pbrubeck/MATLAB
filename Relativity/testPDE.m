function [] = testPDE(m, n)
% Solves a second-order, self-adjoint differential eqation
if nargin<2
    n=m;
end

% Exact solution
f=@(z) real(sin(z.^4));
fx=@(z) 0*real(z);
fy=@(z) 0*real(z);


% Differential operators
E1=eye(m);
E2=eye(n);
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);


% Boundary condtions
a=[1,1;1,1];
b=[0,0;0,0];
B1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
B2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);


% Equation coefficients
A11=1+0*xx.^2.*yy.^2;
A12=0*ones(m,n);
A21=0*ones(m,n);
A22=1+0*xx.^2.*yy.^2;
A3=0*pi^2*sin(pi*(5*xx+3*yy));


% Differential operator
opA=@(uu) Dx*(A11.*(Dx*uu)+A12.*(uu*Dy'))+(A21.*(Dx*uu)+A22.*(uu*Dy'))*Dy'+A3.*uu;


% Right hand sides
xb=[x+1i, x-1i];
yb=[1+1i*y; -1+1i*y];
b1=diag(a(1,:))*f(yb)+diag(b(1,:))*fx(yb);
b2=f(xb)*diag(a(2,:))+fy(xb)*diag(b(2,:));
F=opA(f(xx+1i*yy));


% Preconditioner
[PX1,PY1]=ptop(Dx,E2,A11,Dx,E2);
[PX2,PY2]=ptop(E1,Dy,A12,Dx,E2);
[PX3,PY3]=ptop(Dx,E2,A21,E1,Dy);
[PX4,PY4]=ptop(E1,Dy,A22,E1,Dy);
[PX5,PY5]=ptop(E1,E2,A3 ,E1,E2);
P1=PX1+PX2+PX3+PX4+PX5;
P2=PY1+PY2+PY3+PY4+PY5;


% Solver
[gf,ps,kd,gb,dL]=elliptic(P1,P2,B1,B2,[1,m],[1,n]);
dL(0,Dx,E2,A11,Dx,E2);
dL(1,E1,Dy,A12,Dx,E2);
dL(1,Dx,E2,A21,E1,Dy);
dL(1,E1,Dy,A22,E1,Dy);
dL(1,E1,E2,A3 ,E1,E2);

afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));
tol=2e-15;
maxit=50;

[uu,~,res,its]=bicgstab(afun,rhs,tol,maxit,pfun,[],uu);
uu=gb(uu)+ub;

display(its);
display(res);

% Interpolation
xq=linspace(-1,1,m);
yq=linspace(-1,1,n);
[xxx,yyy]=ndgrid(xq, yq);
uuu=interpcheb(interpcheb(uu,xq,1),yq,2);

err=norm(uu-f(xx+1i*yy),'inf');
errint=norm(uuu-f(xxx+1i*yyy),'inf');
display(err);
display(errint);

% Plot
figure(1);
surf(xxx,yyy,uuu);
colormap(jet(256));
shading interp;
camlight;
end