function [] = brillwaves(n)
kd=2:n-1;
rd=[1,n];
[D,x]=chebD(n); y=x';
L=20;
x=L*x;
D=D/L;

% Boundary conditions
a=[1,1;1,1];
b=0*[L,-L;L,-L];
b1=[0*y+1; 0*y+1];
b2=[0*x+1, 0*x+1];

D2=D*D;
A1=D2+diag(1./x)*D;
A2=D2;

% Imposition of boundary conditions
E=eye(n);
B1=diag(a(1,:))*E(rd,:)+diag(b(1,:))*D(rd,:);
B2=diag(a(2,:))*E(rd,:)+diag(b(2,:))*D(rd,:);
G1=-B1(:,rd)\B1(:,kd);
G2=-B2(:,rd)\B2(:,kd);

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

% Poincare-Steklov operator
N1=zeros(n,2); N1(rd,:)=inv(B1(:,rd));
N2=zeros(n,2); N2(rd,:)=inv(B2(:,rd));
P1=eye(n)-B1'/(B1*B1')*B1;
P2=eye(n)-B2'/(B2*B2')*B2;
F=(B1*B1')*(b1*B2')+(B1*b2)*(B2*B2');
u0=N1*sylvester(B1*B1', B2*B2', F)*N2';
ub=u0+(N1*b1-u0)*P2'+P1*(b2*N2'-u0);

% Problem data
[rr,zz]=ndgrid(x);
k1=4;
k2=3;
q=rr.^2/k1^2.*exp(-rr.^2/k1^2-zz.^2/k2^2);
Q=-1/8*(D2*q+q*D2');
uu=ub;

figure(1);
h=imagesc(x,x,uu);
colormap(jet(256));
colorbar();
axis square;
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('$\rho$');
ylabel('$z$');
drawnow;

% Gauss Seidel
its=8;
for i=1:its
    eqn=A1*uu+uu*A2'-Q.*uu;
    uu=uu-S1*((W1*eqn(kd,kd)*W2')./LL)*S2';
    set(h,'CData',uu);
    drawnow;
    disp(norm(eqn(kd,kd),'fro')/(n-2));
end
end