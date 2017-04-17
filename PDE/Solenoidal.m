function [] = Solenoidal(n)
% Imposes solenoidal and normal boundary conditions to a vector field

% Differential operators
rd=[1,n];
kd=2:n-1;
[D,x]=chebD(n); D2=D*D;
[xx,yy]=ndgrid(x); y=x';

% Original field
zz=xx+1i*yy;
fz=(zz-1/2).*(zz+1/2);
v1=real(fz);
v2=imag(fz);
v1old=v1;
v2old=v2;
vvold=hypot(v1,v2);

% Boundary conditions
a=[0,0;0,0];
b=[1,1;1,1];
b1=v1(rd,:);
b2=v2(:,rd);

% Imposition of boundary conditions
E=eye(n);
B1=diag(a(1,:))*E(rd,:)+diag(b(1,:))*D(rd,:);
B2=diag(a(2,:))*E(rd,:)+diag(b(2,:))*D(rd,:);
G1=-B1(:,rd)\B1(:,kd);
G2=-B2(:,rd)\B2(:,kd);
A1=D2(kd,kd)+D2(kd,rd)*G1;
A2=D2(kd,kd)+D2(kd,rd)*G2;

V1=zeros(n,n-2); V1(rd,:)=G1; V1(kd,:)=eye(n-2);
V2=zeros(n,n-2); V2(rd,:)=G2; V2(kd,:)=eye(n-2);
N1=zeros(n,2); N1(rd,:)=inv(B1(:,rd));
N2=zeros(n,2); N2(rd,:)=inv(B2(:,rd));

P1=V1/(V1'*V1)*V1'; Q1=eye(n)-P1;
P2=V2/(V2'*V2)*V2'; Q2=eye(n)-P2;
F=b1*Q2'*N2+N1'*Q1*b2;
u0=N1*sylvester(N1'*Q1*N1, N2'*Q2*N2, F)*N2';
ub=u0+(N1*b1-u0)*P2'+P1*(b2*N2'-u0);

% Solution
div=D*v1+v2*D';
rhs=div-D2*ub-ub*D2';
[S1,L1]=eig(A1,'vector');
[S2,L2]=eig(A2,'vector');
[L1,L2]=ndgrid(L1,L2); LL=L1+L2;
LL(abs(LL)<1e-5)=inf;
sln=S1*((S1\rhs(kd,kd)/S2')./LL)*S2';
uu=ub+V1*sln*V2';

v1=v1-D*uu;
v2=v2-uu*D';
vv=hypot(v1,v2);

figure(1); clf;
hold on;
imagesc(x,y,vvold);
colormap(jet(256));
quiver(xx,yy,v1old./vvold,v2old./vvold,'k');
hold off;
axis square manual; 
xlabel('x'); ylabel('y');
xlim([-1,1]); ylim([-1,1]);
drawnow;

figure(2); clf;
hold on;
imagesc(x,y,vv);
colormap(jet(256));
quiver(xx,yy,v1./vv,v2./vv,'k');
hold off;
axis square manual; 
xlabel('x'); ylabel('y');
xlim([-1,1]); ylim([-1,1]);
drawnow;

figure(3);
surf(xx,yy,uu);
colormap(jet(256)); 
camlight; shading interp;
axis square manual; 
xlabel('x'); ylabel('y');
drawnow;
end