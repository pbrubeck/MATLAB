function [] = PoissonBinary(m, n, R1, R2, L)
% Solves Poisson for a piecewise, smooth, axisymmetric source defined on 
% the interior and exterior of two spheres 

%% Domain decomposition
% Omega_1: Sphere radius R1, center (0, 0, +L), spherical coords
% Omega_2: Sphere radius R2, center (0, 0, -L), spherical coords
% Omega_3: Complement of Omega_1 + Omega_2, bispherical coords

% Bispherical parameters
a0=sqrt((R1^2-R2^2)^2-8*L^2*(R1^2+R2^2)+16*L^4)/(4*L);
xi1=-a0/R2;
xi2= a0/R1;

m([1,2,3])=m; m1=m(1); m2=m(2); m3=m(3);
E1=eye(m1);
E2=eye(m2);
E3=eye(m3);
E4=eye(n);

% Spherical domains, radial symmetry (0, 1] -> (0, R_i]
% x_i=r/R_i

% Omega_1
[Dx1,x1]=chebD(2*m1);
A1=(Dx1+diag(2./x1))*Dx1;
[A1,Dx1,x1]=radial(A1,Dx1,x1);
g1=1; rd1=1; kd1=2:m1;
a1=1; b1=0;

% Omega_2
[Dx2,x2]=chebD(2*m2);
A2=(Dx2+diag(2./x2))*Dx2;
[A2,Dx2,x2]=radial(A2,Dx2,x2);
g2=1; rd2=1; kd2=2:m2;
a2=1; b2=0;

% Bispherical domain [-1,1] -> [xi1, xi2]
% x3=sinh(xi)

% Omega_3
[Dx3,x3]=chebD(m3);
x3=(xi2-xi1)/2*(x3+1)+xi1; 
Dx3=2/(xi2-xi1)*Dx3;
A3=(diag(x3.^2+1)*Dx3+diag(x3))*Dx3-E3/4;
g3=[1,m3]; rd3=[1,m3]; kd3=1:m3-1;
a3=[1,1]; b3=[0,0];

% Same polar coordinate y=cos(theta) in [-1,1]
[Dy,y]=chebD(n); y=y';
B=(diag(1-y.^2)*Dy+diag(-2*y))*Dy;
rd4=[1,n]; kd4=2:n-1;
a4=[0,0];  b4=[1,1];

%% Coordinate Transformation

% Mobius transformation
M=@(y,a,b) (a*y+b)./(b*y+a);
y1= M(y, R2^2-R1^2-4*L^2, 4*L*R1);
y2=-M(y, R1^2-R2^2-4*L^2, 4*L*R2);

% Tensor product grid
[xx1,yy1]=ndgrid(x1,y);
[xx2,yy2]=ndgrid(x2,y);
[xx3,yy3]=ndgrid(x3,y);

% Physical coordinates (Local to each domain)
r1=R1*x1;
r2=R2*x2;

% Cartesian grids (Midpoint origin)
zzz1=r1*y+L;
xxx1=r1*sqrt(1-y.^2);
zzz2=r2*y-L;
xxx2=r2*sqrt(1-y.^2);
zzz3=a0*(xx3./(sqrt(xx3.^2+1)-yy3))-(sqrt(a0^2+R1^2)-L);
xxx3=a0*sqrt(abs(1-yy3.^2))./(sqrt(xx3.^2+1)-yy3);

%% Equations to solve

% Differential operators
OP1=@(uu) A1*uu+uu*B';
OP2=@(uu) A2*uu+uu*B';
OP3=@(uu) A3*uu+uu*B';

% Interior degrees of freedom
K1=@(uu) uu(kd1,kd4);
K2=@(uu) uu(kd2,kd4);
K3=@(uu) uu(kd3,kd4);

% Right hand side
F1=sin(pi*xx1)./(pi*xx1);
F2=sin(pi*xx2)./(pi*xx2);
F3=zeros(m3,n);

F3=a0^2*((sqrt(xx3.^2+1)-yy3).^(-5/2)).*F3;


%% Boundary Conditions

% Constraint opertor
C1=diag(a1)*E1(rd1,:)+diag(b1)*Dx1(rd1,:);
C2=diag(a2)*E2(rd2,:)+diag(b2)*Dx2(rd2,:);
C3=diag(a3)*E3(rd3,:)+diag(b3)*Dx3(rd3,:);
C4=diag(a4)*E4(rd4,:)+diag(b4)*Dy(rd4,:);

%% Eigenfunctions
function [L,S,G]=eigenfunctions(A,C,m,rd,kd)
    % Give-back matrix
    G=-C(:,rd)\C(:,kd);
    % Eigenfunctions
    S=zeros(m,length(kd));
    [S(kd,:),L]=eig(A(kd,kd)+A(kd,rd)*G,'vector');
    S(rd,:)=G*S(kd,:);
end

[L1,S1,G1]=eigenfunctions(A1,C1,m1,rd1,kd1);
[L2,S2,G2]=eigenfunctions(A2,C2,m2,rd2,kd2);
[L3,S3,G3]=eigenfunctions(A3,C3,m3,rd3,kd3);
[L4,S4,G4]=eigenfunctions(B ,C4,n ,rd4,kd4);

% Interpolated eigenfunctions at the interface
Q1=interpcheb(S4,y1,1);
Q2=interpcheb(S4,y2,1);
Q1(rd4,:)=G4*Q1(kd4,:);
Q2(rd4,:)=G4*Q2(kd4,:);

%% Interior Green functions

[Lx,Ly]=ndgrid(L1,L4); LL1=1./(Lx+Ly);
[Lx,Ly]=ndgrid(L2,L4); LL2=1./(Lx+Ly);
[Lx,Ly]=ndgrid(L3,L4); LL3=1./(Lx+Ly);
GF1=@(rhs) S1*(LL1.*(S1(kd1,:)\rhs/S4(kd4,:).'))*S4.';
GF2=@(rhs) S2*(LL2.*(S2(kd2,:)\rhs/S4(kd4,:).'))*S4.';
GF3=@(rhs) S3*(LL3.*(S3(kd3,:)\rhs/S4(kd4,:).'))*S4.';

% TODO: Account for separation factor
% TODO: Include interpolation

%% Schur complement method

t1=-(Dx1(g1,kd1)*S1(kd1,:)).*((S1(kd1,:)\A1(kd1,g1)).');
t2= (Dx2(g2,kd2)*S2(kd2,:)).*((S2(kd2,:)\A2(kd2,g2)).');
t3=-(Dx3(g3,kd3)*S3(kd3,:)).*((S3(kd3,:)\A3(kd3,g3)).');
t4=-(Dx3(g3,kd3)*S3(kd3,:)).*((S3(kd3,:)\A3(kd3,fliplr(g3))).');

Sig11=S4(kd4,:)*diag(-Dx1(g1,g1)-Dx3(g3(1),g3(1))-t1*LL1-t3(1,:)*LL3)/S4(kd4,:);
Sig22=S4(kd4,:)*diag( Dx2(g2,g2)-Dx3(g3(2),g3(2))-t2*LL2-t3(2,:)*LL3)/S4(kd4,:);
Sig12=S4(kd4,:)*diag(-Dx3(g3(1),g3(2))-t4(1,:)*LL3)/S4(kd4,:);
Sig21=S4(kd4,:)*diag(-Dx3(g3(2),g3(1))-t4(2,:)*LL3)/S4(kd4,:);
Sig=[Sig11, Sig12; Sig21, Sig22];

f1= Dx1(g1,kd1)*K1(GF1(K1(F1)))+Dx3(g3(1),kd3)*K3(GF3(K3(F3)));
f2=-Dx2(g2,kd2)*K2(GF2(K2(F2)))+Dx3(g3(2),kd3)*K3(GF3(K3(F3)));
f=[f1,f2];

ub=zeros(2,n);
ub(:,kd4)=reshape(Sig\f(:),[],2)';
ub(:,rd4)=ub(:,kd4)*G4';

%% Subdomain solutions
uu1=GF1(K1(F1-OP1(E1(:,g1)*ub(1,:))))+E1(:,g1)*ub(1,:);
uu2=GF2(K2(F2-OP2(E2(:,g2)*ub(2,:))))+E2(:,g2)*ub(2,:);
uu3=GF3(K3(F3-OP3(E3(:,g3)*ub)))+E3(:,g3)*ub;
uuu=[uu1;uu2;uu3];

%% Plot Solution
figure(1); clf; hold on;

% Plot Omega_1
surf(xxx1,zzz1,uu1);

% Plot Omega_2
surf(xxx2,zzz2,uu2);

% Grid Clamping
d=2*L+R1+R2;
i1=zzz3>=2*d;
i2=zzz3<=-2*d;
i3=xxx3>=4*d;
[~,imin]=min(abs(x3));
[~,jmin]=max(1-sum(i1|i2|i3));

% Plot Omega_3
ix=1:size(xx3,1); iy=jmin:size(yy3,2);
surf(xxx3(ix,iy),zzz3(ix,iy),uu3(ix,iy));
ix=1:imin-1; iy=1:jmin;
surf(xxx3(ix,iy),zzz3(ix,iy),uu3(ix,iy));
ix=imin+1:size(xx3,1); iy=1:jmin;
surf(xxx3(ix,iy),zzz3(ix,iy),uu3(ix,iy));
hold off;

axis manual;
xlim([0,2*d]); ylim([-d, d]);
colormap(jet(256));
shading interp; camlight;
xlabel('\rho'); ylabel('z');
daspect([1 1 1/sqrt(2)*(max(uuu(:))-min(uuu(:)))/d]);

end

