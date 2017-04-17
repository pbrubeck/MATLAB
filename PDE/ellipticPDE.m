function [] = ellipticPDE(n)
% Solves Laplace's equation with mixed boundary conditions.
% Same as nullBC but optimized.

% Differential operators
rd=[1,n];
kd=2:n-1;
[D,x]=chebD(n); D2=D*D;
[xx,yy]=ndgrid(x); y=x';

% Boundary conditions
a=[1,0;1,0];
b=[0,1;0,1];
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];

% Imposition of boundary conditions
[A1,G1,H1,B1]=setBC(D2,D,a(1,:),b(1,:));
[A2,G2,H2,B2]=setBC(D2,D,a(2,:),b(2,:));

% Eigenfunctions
V1=zeros(n,n-2);
V2=zeros(n,n-2);
[V1(kd,:),L1]=eig(A1,'vector');
[V2(kd,:),L2]=eig(A2,'vector');
V1(rd,:)=G1*V1(kd,:);
V2(rd,:)=G2*V2(kd,:);

% Eigenvalues
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;

% Nullspace
N1=null(D2);
N2=null(D2);

b1=(B1*N1)\b1;
b2=b2/(B2*N2)';
P1=V1/(V1'*V1)*V1'; Q1=eye(n)-P1;
P2=V2/(V2'*V2)*V2'; Q2=eye(n)-P2;
F=b1*Q2'*N2+N1'*Q1*b2;
u0=N1*sylvester(N1'*Q1*N1, N2'*Q2*N2, F)*N2';
ub=u0+(N1*b1-u0)*P2'+P1*(b2*N2'-u0);

% Solution
rhs=-D2*ub-ub*D2';
uu=ub+V1*((V1(kd,:)\rhs(kd,kd)/V2(kd,:)')./LL)*V2';

figure(1);
surf(xx,yy,uu);
colormap(jet(256)); 
camlight; shading interp;
axis square manual; 
xlabel('x'); ylabel('y');
end