function [] = bh(n)
% Black-hole

% Simulation parameters
a0=1;
p=2;

% Differentiation matrix
kd=2:n-1;
rd=[1,n];
[D,x]=chebD(n); y=x';

% Boundary conditions
a=[-1/4,1;0,0];
b=[1,0;1,1];
b1=[0*y; 0*y+1];
b2=[0*x, 0*x];

A1=diag((x+1).^2)*D*D;
A2=diag(1-y.^2)*D*D-diag(2*y)*D;
R=3/256*(p/a0)^2*diag((x+1).^2.*(3-2*x-x.^2).^2);

% Imposition of boundary conditions
E=eye(n);
B1=diag(a(1,:))*E(rd,:)+diag(b(1,:))*D(rd,:);
B2=diag(a(2,:))*E(rd,:)+diag(b(2,:))*D(rd,:);
G1=-B1(:,rd)\B1(:,kd);
G2=-B2(:,rd)\B2(:,kd);

% Poincare-Steklov operator
N1=zeros(n,2); N1(rd,:)=inv(B1(:,rd));
N2=zeros(n,2); N2(rd,:)=inv(B2(:,rd));
P1=eye(n)-B1'/(B1*B1')*B1;
P2=eye(n)-B2'/(B2*B2')*B2;
F=(B1*B1')*(b1*B2')+(B1*b2)*(B2*B2');
u0=N1*sylvester(B1*B1', B2*B2', F)*N2';
ub=u0+(N1*b1-u0)*P2'+P1*(b2*N2'-u0);

% Eigenfunctions
S1=zeros(n,n-2);
S2=zeros(n,n-2);
[S1(kd,:),L1]=eig(A1(kd,kd)+A1(kd,rd)*G1,'vector');
[S2(kd,:),L2]=eig(A2(kd,kd)+A2(kd,rd)*G2,'vector');
S1(rd,:)=G1*S1(kd,:);
S2(rd,:)=G2*S2(kd,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2; LL(abs(LL)<1e-9)=inf;
W1=inv(S1(kd,:));
W2=inv(S2(kd,:));

eqn=A1*ub+ub*A2';
uu=ub-S1*((W1*eqn(kd,kd)*W2')./LL)*S2';
% Gauss Seidel
its=20;
for i=1:its
    eqn=A1*uu+uu*A2'+R*(uu.^(-7));
    uu=uu-S1*((W1*eqn(kd,kd)*W2')./LL)*S2';
    disp(norm(eqn(kd,kd),'fro')/(n-2));
end

% Coordinate mapping
r=(2*a0)./(x+1);
th=acos(x);
rho=r*sin(th)';
z=r*cos(th)';
L=10;

figure(1);
surf(kron([-1,1],rho),z(:,[n:-1:1,1:n]),uu(:,[n:-1:1,1:n]));
xlim([-L,L]);
ylim([-L,L]);
colormap(jet(256));
camlight;
shading interp;
axis square;

E=sqrt(p^2+4*a0^2);
figure(2);
v=(1+2*E./r+6*(a0./r).^2+2*a0^2*E./r.^3+(a0./r).^4).^(1/4);
plot(x, uu(:,1)-v, 'Linewidth', 1);
end