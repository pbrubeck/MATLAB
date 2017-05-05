function [] = PoincareSteklov(n)

% Differential operators
rd=[1,n];
kd=2:n-1;
[D,x]=chebD(n); D2=D*D;
[xx,yy]=ndgrid(x); y=x';

% Boundary conditions
a=[1,1;1,1];
b=[0,0;0,0];
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];

% Imposition of boundary conditions
E=eye(n);
B1=diag(a(1,:))*E(rd,:)+diag(b(1,:))*D(rd,:);
B2=diag(a(2,:))*E(rd,:)+diag(b(2,:))*D(rd,:);
G1=-B1(:,rd)\B1(:,kd);
G2=-B2(:,rd)\B2(:,kd);
A1=D2(kd,kd)+D2(kd,rd)*G1;
A2=D2(kd,kd)+D2(kd,rd)*G2;

% Eigenvectors
S1=zeros(n,n-2);
S2=zeros(n,n-2);
[S1(kd,:),L1]=eig(A1,'vector');
[S2(kd,:),L2]=eig(A2,'vector');
S1(rd,:)=G1*S1(kd,:);
S2(rd,:)=G2*S2(kd,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;
LL(abs(LL)<1e-9)=inf;

% Poincare-Steklov operator
N1=zeros(n,2); N1(rd,:)=inv(B1(:,rd));
N2=zeros(n,2); N2(rd,:)=inv(B2(:,rd));
P1=eye(n)-B1'/(B1*B1')*B1;
P2=eye(n)-B2'/(B2*B2')*B2;
F=(B1*B1')*(b1*B2')+(B1*b2)*(B2*B2');
u0=N1*sylvester(B1*B1', B2*B2', F)*N2';
ub=u0+(N1*b1-u0)*P2'+P1*(b2*N2'-u0);

% Solution
rhs=-D2*ub-ub*D2';
uu=ub+S1*((S1(kd,:)\rhs(kd,kd)/S2(kd,:)')./LL)*S2';

figure(1);
surf(xx,yy,uu);
colormap(jet(256)); 
camlight; shading interp;
axis square manual; 
xlabel('x'); ylabel('y');
end