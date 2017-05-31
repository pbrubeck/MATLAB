function [] = PoincareSteklov(n)
% Differential operators
n([1,2])=n;
rd1=[1,n(1)];
kd1=2:n(1)-1;
rd2=[1,n(2)];
kd2=2:n(2)-1;
[Dx,x]=chebD(n(1)); A1=Dx*Dx;
[Dy,y]=chebD(n(2)); A2=Dy*Dy;
[xx,yy]=ndgrid(x,y); y=y';

% Boundary conditions
a=[1,0;1,0];
b=[0,1;0,1];
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];

% Imposition of boundary conditions
E1=eye(n(1));
E2=eye(n(2));
B1=diag(a(1,:))*E1(rd1,:)+diag(b(1,:))*Dx(rd1,:);
B2=diag(a(2,:))*E2(rd2,:)+diag(b(2,:))*Dy(rd2,:);
G1=-B1(:,rd1)\B1(:,kd1);
G2=-B2(:,rd2)\B2(:,kd2);
S1=A1(kd1,kd1)+A1(kd1,rd1)*G1;
S2=A2(kd2,kd2)+A2(kd2,rd2)*G2;

% Eigenvectors
V1=zeros(n(1),n(1)-2);
V2=zeros(n(2),n(2)-2);
[V1(kd1,:),L1]=eig(S1,'vector');
[V2(kd2,:),L2]=eig(S2,'vector');
V1(rd1,:)=G1*V1(kd1,:);
V2(rd2,:)=G2*V2(kd2,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;
LL(abs(LL)<1e-9)=inf;

% Poincare-Steklov operator
N1=zeros(n(1),2); N1(rd1,:)=inv(B1(:,rd1));
N2=zeros(n(2),2); N2(rd2,:)=inv(B2(:,rd2));
P1=eye(n(1))-B1'/(B1*B1')*B1;
P2=eye(n(2))-B2'/(B2*B2')*B2;
F=(B1*B1')*(b1*B2')+(B1*b2)*(B2*B2');
u0=N1*sylvester(B1*B1', B2*B2', F)*N2';
ub=u0+(N1*b1-u0)*P2+P1*(b2*N2'-u0);

% Solution
rhs=-A1*ub-ub*A2';
uu=ub+V1*((V1(kd1,:)\rhs(kd1,kd2)/V2(kd2,:)')./LL)*V2';

figure(1);
surf(xx,yy,uu);
colormap(jet(256)); 
camlight; shading interp;
axis square manual; 
xlabel('x'); ylabel('y');
end