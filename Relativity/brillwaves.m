function [] = brillwaves(n)
% Black-hole plus brill wave

% Simulation parameters
a0=2;
A0=1;
eta0=0.5;
s0=1;
m=4;
L=24;

% Differentiation matrix
n([1,2])=n;
[Dx,x]=chebD(n(1));
[Dy,y]=chebD(n(2)); y=y';

% Boundary conditions
a=[-1/4,1;0,0];
b=[1,0;1,1];
b1=[0*y; 0*y+1];
b2=[0*x, 0*x];

A1=diag((x+1).^2)*Dx*Dx;
A2=diag(1-y.^2)*Dy*Dy-diag(2*y)*Dy;

% Imposition of boundary conditions
E1=eye(n(1));
E2=eye(n(2));
B1=diag(a(1,:))*E1([1,end],:)+diag(b(1,:))*Dx([1,end],:);
B2=diag(a(2,:))*E2([1,end],:)+diag(b(2,:))*Dy([1,end],:);
G1=-B1(:,[1,end])\B1(:,2:end-1);
G2=-B2(:,[1,end])\B2(:,2:end-1);

% Poincare-Steklov operator
N1=zeros(n(1),2); N1([1,end],:)=inv(B1(:,[1,end]));
N2=zeros(n(2),2); N2([1,end],:)=inv(B2(:,[1,end]));
P1=eye(n(1))-B1'/(B1*B1')*B1;
P2=eye(n(2))-B2'/(B2*B2')*B2;
F=(B1*B1')*(b1*B2')+(B1*b2)*(B2*B2');
u0=N1*sylvester(B1*B1', B2*B2', F)*N2';
ub=u0+(N1*b1-u0)*P2'+P1*(b2*N2'-u0);

% Eigenfunctions
S1=zeros(n(1),n(1)-2);
S2=zeros(n(2),n(2)-2);
[S1(2:end-1,:),L1]=eig(A1(2:end-1,2:end-1)+A1(2:end-1,[1,end])*G1,'vector');
[S2(2:end-1,:),L2]=eig(A2(2:end-1,2:end-1)+A2(2:end-1,[1,end])*G2,'vector');
S1([1,end],:)=G1*S1(2:end-1,:);
S2([1,end],:)=G2*S2(2:end-1,:);
[L1,L2]=ndgrid(L1,L2);
LL=L1+L2; LL(abs(LL)<1e-9)=inf;
W1=inv(S1(2:end-1,:));
W2=inv(S2(2:end-1,:));

eta=log(2./(x+1));
qq=A0*(exp(-((eta+eta0)/s0).^2)+exp(-((eta-eta0)/s0).^2))*((1-y.^2).^(m/2));
A3=diag((x+1).^2)*Dx*Dx+diag(x+1)*Dx;
A4=diag(1-y.^2)*Dy*Dy-diag(y)*Dy;
R=(A3*qq+qq*A4')/4;

eqn=A1*ub+ub*A2';
uu=ub-S1*((W1*eqn(2:end-1,2:end-1)*W2')./LL)*S2';
% Gauss Seidel
its=90;
for i=1:its
    eqn=A1*uu+uu*A2'+R.*uu;
    uu=uu-S1*((W1*eqn(2:end-1,2:end-1)*W2')./LL)*S2';
    disp(norm(eqn(2:end-1,2:end-1),'fro')/sqrt(prod((n-2))));
end

% Coordinate mapping
r=(2*a0)./(x+1);
th=acos(y);
rho=r*sin(th);
z=r*cos(th);

figure(1);
mesh(kron([-1,1],rho),z(:,[n(2):-1:1,1:n(2)]), uu(:,[n(2):-1:1,1:n(2)]));
colormap(jet(256));
%camlight; shading interp;
axis square;
xlim([-L,L]);
ylim([-L,L]);
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$\rho$');
ylabel('$z$');
zlabel('$\psi$');
end