function [] = Poisson2sphere(m,n,R0)
% Solves Poisson for a source that vanishes outside a spherical region
% Domain decomposition: Unit sphere + spherical shell to infinity

% Two radial coordinates r in (0, R0] U [R0, inf)
m([1,2])=m; m1=m(1); m2=m(2); 
E1=eye(m1);
E2=eye(m2);
E3=eye(n);

% Inner sphere, radial symmetry (0, 1] -> (0, R0]
% x1=r/R0
[Dx1,x1]=chebD(2*m1);
A1=(Dx1+diag(2./x1))*Dx1;
[A1,Dx1,x1]=radial(A1,Dx1,x1);
g1=1;
rd1=1;
kd1=2:m1;
a1=1; 
b1=0;

% Outer shell, vanishing at infinity [-1, 1] -> [R0, inf)
% x2=2*R0/r-1
[Dx2,x2]=chebD(m2);
A2=diag((x2+1).^2)*Dx2*Dx2;
g2=1;
rd2=[1,m2];
kd2=2:m2-1;
a2=[1,1];
b2=[0,0];

% Same polar coordinate y=cos(theta) in [-1,1]
[Dy,y]=chebD(n);
y=y';
B=(diag(1-y.^2)*Dy+diag(-2*y))*Dy;
rd3=[1,n];
kd3=2:n-1;
a3=[0,0]; 
b3=[1,-1];
K1=@(uu) uu(kd1,kd3);
K2=@(uu) uu(kd2,kd3);

% Differential operators
OP1=@(uu) A1*uu+uu*B';
OP2=@(uu) A2*uu+uu*B';

% Cooridinate grid
r1=R0*x1;
r2=2*R0./(x2+1);
r=[r2;r1];
[rr1,th1]=ndgrid(r1,acos(y));
[rr2,th2]=ndgrid(r2,acos(y));
rr=[rr2;rr1];
th=[th2;th1];

% Source
F1=sin(pi*rr1/R0)./(pi*rr1/R0);
F2=zeros(m2,n);

% Constraint operator
C1=diag(a1)*E1(rd1,:)+diag(b1)*Dx1(rd1,:);
C2=diag(a2)*E2(rd2,:)+diag(b2)*Dx2(rd2,:);
C3=diag(a3)*E3(rd3,:)+diag(b3)*Dy(rd3,:);

% Giveback matrix
G1=-C1(:,rd1)\C1(:,kd1);
G2=-C2(:,rd2)\C2(:,kd2);
G3=-C3(:,rd3)\C3(:,kd3);

% Interior Schur complement
SA1=A1(kd1,kd1)+A1(kd1,rd1)*G1;
SA2=A2(kd2,kd2)+A2(kd2,rd2)*G2;
SB=B(kd3,kd3)+B(kd3,rd3)*G3;

% Eigendecomposition
S1=zeros(m1,length(kd1));
S2=zeros(m2,length(kd2));
S3=zeros(n ,length(kd3));
[S1(kd1,:),L1]=eig(SA1,'vector');
[S2(kd2,:),L2]=eig(SA2,'vector');
[S3(kd3,:),L3]=eig(SB ,'vector');

S1(rd1,:)=G1*S1(kd1,:);
S2(rd2,:)=G2*S2(kd2,:);
S3(rd3,:)=G3*S3(kd3,:);

% Interior Green functions
[Lx,Ly]=ndgrid(L1,L3); LL1=1./(Lx+Ly);
[Lx,Ly]=ndgrid(L2,L3); LL2=1./(Lx+Ly);
GF1=@(rhs) S1*(LL1.*(S1(kd1,:)\rhs/S3(kd3,:).'))*S3.';
GF2=@(rhs) S2*(LL2.*(S2(kd2,:)\rhs/S3(kd3,:).'))*S3.';

% Interface Schur complement
f1= 1*(S1(kd1,:)\A1(kd1,g1))'.*(Dx1(g1,kd1)*S1(kd1,:));
f2=-2*(S2(kd2,:)\A2(kd2,g2))'.*(Dx2(g2,kd2)*S2(kd2,:));
SS=S3(kd3,:)*diag((Dx1(g1,g1)-2*Dx2(g2,g2)-f1*LL1-f2*LL2))/S3(kd3,:);

% Naive computation of Schur Complement
function v=schurcomp(rhs)
    gb=zeros(1,n);
    gb(kd3)=rhs;
    gb(rd3)=gb(kd3)*G3';
    f1=K1(OP1(E1(:,g1)*gb));
    f2=K2(OP2(E2(:,g2)*gb));
    v=(Dx1(g1,g1)-2*Dx2(g2,g2))*rhs-Dx1(g1,kd1)*K1(GF1(f1))+2*Dx2(g2,kd2)*K2(GF2(f2));
end
if false
SS2=eye(n-2);
for i=1:n-2
    SS2(i,:)=schurcomp(SS2(i,:));
end
SS2=SS2.';
figure(3);
imagesc(SS-SS2); colorbar();
end

% Solve for the interface
rhs=Dx1(g1,kd1)*K1(GF1(K1(F1)))-2*Dx2(g2,kd2)*K2(GF2(K2(F2)));
ug=zeros(1,n);
ug(kd3)=rhs/SS';
ug(rd3)=ug(:,kd3)*G3';

% Solve for each domain
u1=GF1(K1(F1-OP1(E1(:,g1)*ug)))+E1(:,g1)*ug;
u2=GF2(K2(F2-OP2(E2(:,g2)*ug)))+E2(:,g2)*ug;
uu=[u2;u1];

% Exact solution
phi0=@(r) R0^2/pi^2*((r>=R0).*(-R0./(r+(r==0)))+(r<R0).*(-1-sinc(pi/R0*r)));

% Cartesian grid
xx=r*sqrt(1-y.^2);
zz=r*y;

err=uu-phi0(rr);
err(isnan(err))=0;
err=norm(err,'inf');
display(err);

% Plot
figure(1);
surf(xx,zz,uu);
colormap(jet(256));
colorbar;
view(2); 
d=4*R0; xlim([0,d]); ylim([-d,d]);
daspect([1 1 1/sqrt(2)*(max(uu(:))-min(uu(:)))/d]);
shading interp; camlight;
xlabel('\rho'); ylabel('z');
end