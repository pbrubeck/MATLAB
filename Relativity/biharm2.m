function [] = biharm2(m,n)
% Biharmonic equation in two dimensions
if(nargin<2)
    n=m;
end

% Exact solution
f=@(z) real(sin(z.^4));
fx=@(z) real(4*z.^3.*cos(z.^4));
fy=@(z) real(1i*4*z.^3.*cos(z.^4));

% Differential operators
E1=eye(m);
E2=eye(n);
[Dx,x]=chebD(m); Dxx=Dx*Dx;
[Dy,y]=chebD(n); Dyy=Dy*Dy; y=y';
[xx,yy]=ndgrid(x,y);


% Boundary operators
rd1=[1,2,m-1,m]; % Removed dof at x
rd2=[1,2,n-1,n]; % Removed dof at y
be1=[1,1,m,m];   % Boundary eqn at x
be2=[1,1,n,n];   % Boundary eqn at y
a=[1,0,0,1;1,0,0,1];
b=[0,1,1,0;0,1,1,0];
c=[0,0,0,0;0,0,0,0];
B1=diag(a(1,:))*E1(be1,:)+diag(b(1,:))*Dx(be1,:)+diag(c(1,:))*Dxx(be1,:);
B2=diag(a(2,:))*E2(be2,:)+diag(b(2,:))*Dy(be2,:)+diag(c(1,:))*Dyy(be2,:);

% Differential operator
lap=@(uu) Dxx*uu+uu*Dyy';
opA=@(uu) lap(lap(uu));


A11=ones(m,n);
A12=ones(m,n);
A21=ones(m,n);
A22=ones(m,n);

% Right hand sides
zz=xx+1i*yy;
z1=zz(be1,:);
z2=zz(:,be2);
b1=diag(a(1,:))*f(z1)+diag(b(1,:))*fx(z1);
b2=f(z2)*diag(a(2,:))+fy(z2)*diag(b(2,:));
F=opA(f(xx+1i*yy));

% Preconditioner
P1=Dxx*Dxx;
P2=Dyy*Dyy;

% Solver
[gf,ps,kd,gb,dL]=elliptic(P1,P2,B1,B2,rd1,rd2);
dL(0,Dxx,E2,A11,Dxx,E2);
dL(1,E1,Dyy,A12,Dxx,E2);
dL(1,Dxx,E2,A21,E1,Dyy);
dL(1,E1,Dyy,A22,E1,Dyy);


afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));
tol=8*eps;
maxit=100;

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

