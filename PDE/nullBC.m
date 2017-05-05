function [] = nullBC(N)
% Solves Laplace's equation with inhomogeneous boundary conditions

% Differential operators
rd=[1,N];
kd=2:N-1;
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x); y=x';

% Boundary conditions
a=[1,0;1,0];
b=[0,1;0,1];
b1=[0.2*sin(3*pi*y); 0*y];
b2=[(x<0).*sin(pi*x).^4, 0*x];

V1=zeros(N); L1=zeros(N,1);
V2=zeros(N); L2=zeros(N,1);

% Nullspace
V1(:,rd)=null(D2);
V2(:,rd)=null(D2);

% Imposition of boundary conditions
[A1,G1,H1,C1]=setBC(D2,D,a(1,:),b(1,:));
[A2,G2,H2,C2]=setBC(D2,D,a(2,:),b(2,:));

% Eigenfunctions
[V1(kd,kd),L1(kd)]=eig(A1,'vector');
[V2(kd,kd),L2(kd)]=eig(A2,'vector');
V1(rd,kd)=G1*V1(kd,kd);
V2(rd,kd)=G2*V2(kd,kd);

% Eigenvalues
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2;

% Boundary operator on nullspace
P11=kron(V2(:,kd), C1*V1(:,rd));
P12=zeros(2*N, 2*(N-2));
P13=kron(V2(:,rd), C1*V1(:,rd));
P21=zeros(2*N, 2*(N-2));
P22=kron(C2*V2(:,rd), V1(:,kd));
P23=kron(C2*V2(:,rd), V1(:,rd));
P=[P11, P12, P13; P21, P22, P23];

% Inhomogenous boundary term
qq=P\[b1(:); b2(:)];
U=zeros(N);
U(rd,kd)=reshape(qq(1:2*(N-2)), [2,N-2]);
U(kd,rd)=reshape(qq(1+2*(N-2):4*(N-2)), [N-2,2]);
U(rd,rd)=reshape(qq(1+4*(N-2):4*(N-1)), [2,2]);
ub=V1*U*V2';

% Solution
rhs=-D2*ub-ub*D2';
U(kd,kd)=(V1(kd,kd)\rhs(kd,kd)/V2(kd,kd)')./LL(kd,kd);
uu=V1*U*V2';

figure(1);
surf(xx,yy,uu);
colormap(jet(256)); 
camlight; shading interp;
axis square manual; 
xlabel('x'); ylabel('y');
end