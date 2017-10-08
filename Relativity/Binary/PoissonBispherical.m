function [] = PoissonBispherical(m,n,L,R1,R2)
% Solves Poisson's equation in bispherical coordinats with mixed boundary conditions.

% L : source distance
% R1, R2 : source radii

a0=sqrt((R1^2-R2^2)^2-8*L^2*(R1^2+R2^2)+16*L^4)/(4*L);
x1=-a0/R2;
x2= a0/R1; 

% Set coordinate grid [x,y]=[sinh(xi), cos(eta)]
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
x=(x2-x1)/2*(x+1)+x1; Dx=2/(x2-x1)*Dx;
[xx,yy]=ndgrid(x,y);

% Differential operators
E1=eye(m);
E2=eye(n);
A=(diag(x.^2+1)*Dx+diag(x))*Dx-E1/8;
B=(diag(1-y.^2)*Dy-diag(2*y))*Dy-E2/8;

% Right hand side
zc=a0*(xx./(sqrt(xx.^2+1)-yy)-sqrt(1+x2^-2))+L;
xc=a0*sqrt(abs(1-yy.^2))./(sqrt(xx.^2+1)-yy);
r1=hypot(xc,zc-L);
r2=hypot(xc,zc+L);
src0=@(r) (r<=1).*sinc(pi*r);
src=src0(r1/R1)+src0(r2/R2);
rhs=a0^2*((sqrt(xx.^2+1)-yy).^(-5/2)).*src;

% Analytic solution
phi0=@(r) (1/pi^2)*((r>=1).*(-1./(r+(r==0)))+(r<1).*(-1-sinc(pi*r)));
phi=R1^2*phi0(r1/R1)+R2^2*phi0(r2/R2);

% Kept degrees of freedom
rd1=[1,m]; kd1=2:m-1;
rd2=[1,n]; kd2=2:n-1;
kd=@(uu) uu(kd1,kd2);

% Boundary conditions
a=[1,1;0,0];
b=[0,0;1,1];

% Constraint opertor
C1=diag(a(1,:))*E1(rd1,:)+diag(b(1,:))*Dx(rd1,:);
C2=diag(a(2,:))*E2(rd2,:)+diag(b(2,:))*Dy(rd2,:);

% Boundary data
phi0=phi./sqrt(sqrt(xx.^2+1)-yy);
b1=C1*phi0;
b2=phi0*C2';

% Poincare-Steklov operator
N1=zeros(m,length(rd1)); N1(rd1,:)=inv(C1(:,rd1));
N2=zeros(n,length(rd2)); N2(rd2,:)=inv(C2(:,rd2));
P1=eye(m)-C1'/(C1*C1')*C1;
P2=eye(n)-C2'/(C2*C2')*C2;
u0=N1*sylvester(C1*C1', C2*C2', (C1*C1')*(b1*C2')+(C1*b2)*(C2*C2'))*N2';
ub=u0+(N1*b1-u0)*P2+P1*(b2*N2'-u0);

% Give-back matrix
G1=-C1(:,rd1)\C1(:,kd1);
G2=-C2(:,rd2)\C2(:,kd2);

% Schur complement
SA=A(kd1,kd1)+A(kd1,rd1)*G1;
SB=B(kd2,kd2)+B(kd2,rd2)*G2;

% Eigenfunctions
V1=zeros(m,length(kd1));
V2=zeros(n,length(kd2));
[V1(kd1,:),L1]=eig(SA,'vector');
[V2(kd2,:),L2]=eig(SB,'vector');
V1(rd1,:)=G1*V1(kd1,:);
V2(rd2,:)=G2*V2(kd2,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;

% Solve
uu=ub+V1*(((V1(kd1,:)\kd(rhs-(A*ub+ub*B.'))/V2(kd2,:).'))./LL)*V2.';
uu=sqrt(sqrt(xx.^2+1)-yy).*uu;

% Interpolation
if false
    M=1024;
    N=1024;
    xq=(x2-x1)/2*(cos(pi/(M-1)*(0:M-1))+1)+x1;
    yq=cos(pi/(N-1)*(0:N-1));
    ran=@(x) 2*(x-min(x))/(max(x)-min(x))-1;
    uuu=interpcheb(interpcheb(uu,ran(xq),1),ran(yq),2);
    [xq,yq]=ndgrid(xq,yq);
else
    xq=xx; yq=yy; uuu=uu;
end

% Transformation
z0=sqrt(a0^2+R1^2)-L;
zzz=a0*(xq./(sqrt(xq.^2+1)-yq))-z0;
xxx=a0*sqrt(abs(1-yq.^2))./(sqrt(xq.^2+1)-yq);

% Grid Clamping
d=2*L+R1+R2;
i1=zzz>=2*d;
i2=zzz<=-2*d;
i3=xxx>=4*d;
[~,imin]=min(abs(xq));
[~,jmin]=max(1-sum(i1|i2|i3));

% Plot
figure(1); clf; hold on;
ix=1:size(xq,1); iy=jmin:size(yq,2);
surf(xxx(ix,iy),zzz(ix,iy),uuu(ix,iy));
ix=1:imin-1; iy=1:jmin;
surf(xxx(ix,iy),zzz(ix,iy),uuu(ix,iy));
ix=imin+1:size(xq,1); iy=1:jmin;
surf(xxx(ix,iy),zzz(ix,iy),uuu(ix,iy));

axis manual;
xlim([0,d]); ylim([-d, d]);
daspect([1 1 1/sqrt(2)*(max(uuu(:))-min(uuu(:)))/d]);
colormap(jet(256));
shading interp; camlight;
xlabel('\rho'); ylabel('z');
hold off;
end