function [] = brill(m,n)
if nargin<2
    n=m;
end

% Simulation parameters
L=8;
r0=sqrt(2)*L;
A0=1;

[Dx,x]=chebD(2*m);
A1=(diag(x.^2)*Dx+diag(2*x))*Dx;
[A1,Dx,x]=radial(A1,Dx,x);

[Dy,y]=chebD(n); y=y';
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;

a=[1,1;0,0];
b=[1,1;1,1];

% Imposition of boundary conditions
E1=eye(m);
E2=eye(n);
B1=a(1,1)*E1(1,:)+b(1,1)*Dx(1,:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);

b1=0*y+1;
b2=[0*x,0*x];

% Coordinate mapping
r=r0*x;
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

qq=A0*(rho.^2).*exp(-(rho/1).^2-(z/1).^2);
A3=1/4*((diag(x.^2)*Dx*Dx+diag(x)*Dx)*qq+qq*(diag(1-y.^2)*Dy*Dy-diag(y)*Dy)');
F=zeros(m,n);
opA=@(uu) A1*uu+uu*A2'+A3.*uu;

% Solution
[gf,ps,kd,gb]=elliptic(A1,E2,E1,A2,B1,B2,1,[1,n]);

ub=ps(b1,b2);
rhs=kd(F-opA(ub));
uu=kd(ub+gf(rhs));

afun=@(uu) kd(opA(gb(uu)));
pfun=@(uu) kd(gf(uu));
tol=2e-15;
maxit=50;

[uu,~,res,its]=bicgstab(afun,rhs,tol,maxit,pfun,[],uu);
uu=gb(uu)+ub;

display(its);
display(res);

[rr,zz]=ndgrid(linspace(-L,L,2*m));
xx=hypot(rr,zz)/r0;
yy=zz./hypot(rr,zz);
uuu=interpcheb2(uu([1:end,end:-1:1],:),xx,yy);


figure(1);
surf([rho; -flipud(rho)], z([1:end,end:-1:1],:), uu([1:end,end:-1:1],:));
surf(rr,zz,uuu);

colormap(jet(256));
colorbar;
shading interp; 
%camlight; 
axis square;
xlim([-L,L]);
ylim([-L,L]);
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$\rho$');
ylabel('$z$');
title('$\psi(\rho,z)$');
view(2);
end