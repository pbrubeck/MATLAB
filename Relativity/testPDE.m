function [] = testPDE(m, n)
% Solves a second-order, self-adjoint differential eqation
if nargin<2
    n=m;
end

% Exact solution
f=@(z) real(sin(z.^4));
fx=@(z) real(4*z.^3.*cos(z.^4));
fy=@(z) real(1i*4*z.^3.*cos(z.^4));

% Differential operators
E1=eye(m);
E2=eye(n);
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);


% Boundary condtions
a=[1,1;1,1];
b=[1,1;1,1];
B1=diag(a(1,:))*E1([1,m],:)+diag(b(1,:))*Dx([1,m],:);
B2=diag(a(2,:))*E2([1,n],:)+diag(b(2,:))*Dy([1,n],:);


% Equation coefficients
A11=ones(m,n);
A21=ones(m,n)/3;
A12=ones(m,n)/3;
A22=ones(m,n);
A0=pi^2*sin(pi*(5*xx+3*yy));

if(any((A12+A21).^2>=4*A11.*A22))
    % This costed a day of work
    error('Equation is not elliptic');
end

% Differential operator
ell=@(uu,ux,uy) Dx*(A11.*ux+A12.*uy)+(A21.*ux+A22.*uy)*Dy'+A0.*uu;
opA=@(uu) ell(uu, Dx*uu, uu*Dy');


% Right hand sides
xb=[x+1i, x-1i];
yb=[1+1i*y; -1+1i*y];
b1=diag(a(1,:))*f(yb)+diag(b(1,:))*fx(yb);
b2=f(xb)*diag(a(2,:))+fy(xb)*diag(b(2,:));
F=opA(f(xx+1i*yy));


% Preconditioner: Partial traces of the full operator
[PX1, ~ ]=ptop(Dx,E2,A11,Dx,E2);
[ ~ ,PY2]=ptop(E1,Dy,A22,E1,Dy);
[PX3,PY3]=ptop(E1,E2,A0, E1,E2);
P1=(PX1+PX3)/n; 
P2=(PY2+PY3)/m;
P1=P1-trace(P1)/(2*m)*eye(m);
P2=P2-trace(P2)/(2*n)*eye(n);

% Solver
[gf,ps,kd,gb,dL]=elliptic(P1,P2,B1,B2,[1,m],[1,n]);
dL(0,Dx,E2,A11,Dx,E2);
dL(1,E1,Dy,A12,Dx,E2);
dL(1,Dx,E2,A21,E1,Dy);
dL(1,E1,Dy,A22,E1,Dy);
dL(1,E1,E2,A0, E1,E2);

afun=@(uu)  kd(opA(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));
tol=2*eps;
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