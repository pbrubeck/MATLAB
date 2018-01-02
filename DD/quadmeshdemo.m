function [ ] = quadmeshdemo( m )
% Attempting to combine Schur complement and NKP with quadrileteral
% domains.

% Differential operators
n=m;
[Dx,x]=chebD(m);
[Dy,y]=chebD(n);
% (xx,yy) fine grid for over-integration
[xx,wx]=gaulob(-1,1,m+32);
[yy,wy]=gaulob(-1,1,n+32);
% (xxx,yyy) is two dimensional, here we evaluate the equation coefficients
[xxx,yyy]=ndgrid(xx,yy);

% Vertices
v0=[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)];
%v0=2/(2-sqrt(2))*[1i;0;1];

% Sides
L=abs(v0([3,1,2])-v0([2,3,1]));
% Contact triangle
V=eye(3)+diag((sum(L)/2-L)./L([2,3,1]))*[-1,0,1; 1,-1,0; 0,1,-1];

z0=zeros(7,1);
z0([1,2,3])=v0;       % Vertices
z0([4,5,6])=V*v0;     % Touch points
z0(7)=(L'*v0)/sum(L); % Incenter

Z1=reshape(z0([7,4,5,1]),[2,2]);
Z2=reshape(z0([7,5,6,2]),[2,2]);
Z3=reshape(z0([7,6,4,3]),[2,2]);

% Evaluate jacobian determinant and metric tensor [E, F; F, G]
[zz1,J1,E1,F1,G1]=mapquad(Z1,xxx,yyy);
[zz2,J2,E2,F2,G2]=mapquad(Z2,xxx,yyy);
[zz3,J3,E3,F3,G3]=mapquad(Z3,xxx,yyy);

figure(1);
mesh(real(zz1),imag(zz1),F1); hold on;
mesh(real(zz2),imag(zz2),F2);
mesh(real(zz3),imag(zz3),F3);

hold off;
axis equal;
view(2);

% Galerkin stiffness and mass (matrix-free) operators, with their NKP
[lap1,mass1,Mx1,My1,Kx1,Ky1]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J1,E1,F1,G1);
[lap2,mass2,Mx2,My2,Kx2,Ky2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J2,E2,F2,G2);
[lap3,mass3,Mx3,My3,Kx3,Ky3]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J3,E3,F3,G3);

% Helper routine to solve generalized Sylvester equations
I1=eye(m);
I2=eye(n);
B1=I1([1,m],:);
B2=I2([1,n],:);
[gf,ps,kd,gb]=elliptic(Kx1,My1,Mx1,Ky1,B1,B2,[1,m],[1,n]);
b1=zeros(2,n);
b2=zeros(m,2);
F=ones(m,n);

afun=@(uu)  kd(lap1(gb(uu)));
pfun=@(rhs) kd(gf(rhs));
ub=ps(b1,b2);
rhs=kd(mass1(F)-lap1(ub));
tol=1e-11;
maxit=50;

% Solve the big problem in one domain, the Schur complement is missing
[uu,~,res,its]=bicgstab(afun,rhs,tol,maxit,pfun,[],kd(ub+gf(rhs)));
uu=gb(uu)+ub;

display(its);
display(res);

[xx,yy]=ndgrid(x);
zz0=mapquad(Z1,xx,yy);

figure(3);
surf(real(zz0), imag(zz0), reshape(uu,size(zz0)));
colormap(jet(256));
shading interp; camlight; view(2);
dx=diff(xlim());
dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
end