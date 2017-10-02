function [] = PoissonProlate(m,n,f,zmax,R0)
% Solves Poisson's equation in oblate coordinats with mixed boundary conditions.

% f : semi-focal length
% zmax : semi-minor axis
% R0 : source radius

x0=zmax/f;

% Set coordinate grid [x,y]=[cosh(xi), cos(eta)]
[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
x=(x0-1)/2*(x+1)+1; Dx=2/(x0-1)*Dx;
[xx,yy]=ndgrid(x,y);

% Differential operators
A=(diag(x.^2-1)*Dx+diag(2*x))*Dx;
B=(diag(1-y.^2)*Dy-diag(2*y))*Dy;
E1=eye(m);
E2=eye(n);

% Right hand side
r1=f*(xx-yy);
r2=f*(xx+yy);

src0=@(r) (r<=R0).*sinc(pi/R0*r);
src1=@(r) exp(-r.^2/(2*R0^2));
src2=@(r) (r.^2-3*R0^2).*exp(-r.^2/(2*R0^2))/R0^4;
src=src1(r1)+src1(r2);
rhs=diag(f^2*x.^2)*src-src*diag(f^2*y.^2);

% Analytic solution
phi0=@(r) R0^2/pi^2*((r>=R0).*(-R0./(r+(r==0)))+(r<R0).*(-1-sinc(pi/R0*r)));
phi1=@(r) -R0^3*(sqrt(pi/2)*erf(r/(sqrt(2)*R0))+(r==0)/R0)./(r+(r==0));
phi2=@(r) exp(-r.^2/(2*R0^2));
phi=phi1(r1)+phi1(r2);

% Kept degrees of freedom
rd1=[1,m]; kd1=2:m-1;
rd2=[1,n]; kd2=2:n-1;
kd=@(uu) uu(kd1,kd2);

% Boundary conditions
a=[1,0;0,0];
b=[x0-1,1;1,1];

% Constraint opertor
C1=diag(a(1,:))*E1(rd1,:)+diag(b(1,:))*Dx(rd1,:);
C2=diag(a(2,:))*E2(rd2,:)+diag(b(2,:))*Dy(rd2,:);

% Boundary data
b1=C1*phi;
b2=phi*C2';

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
% Interpolation
if false
    M=1024;
    N=1024;
    xq=(x0-1)/2*(cos(pi/(M-1)*(0:M-1))+1)+1;
    yq=cos(pi/(N-1)*(0:N-1));
    ran=@(x) 2*(x-min(x))/(max(x)-min(x))-1;
    uuu=interpcheb(interpcheb(uu,ran(xq),1),ran(yq),2);
    [xq,yq]=ndgrid(xq,yq);
else
    xq=xx; yq=yy; uuu=uu;
end

% Transformation
zzz=f*xq.*yq;
xxx=f*sqrt((xq.^2-1).*(1-yq.^2));
figure(1);
surf([xxx;-flipud(xxx)],zzz([1:end, end:-1:1],:),uuu([1:end, end:-1:1],:));
colormap(jet(256)); 
camlight; shading interp; 
view(2); axis manual;
daspect([1 1 sqrt(2)*(max(uuu(:))-min(uuu(:)))/zmax]);
xlabel('\rho'); ylabel('z');

figure(2);
plot(y,C1*uu);

figure(3);
plot(x,uu*C2');
end