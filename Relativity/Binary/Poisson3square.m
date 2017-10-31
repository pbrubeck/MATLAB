function [] = Poisson3square(m, n)
% Solves Poisson for three connected regions, same topology as
% PoissonBinary

%% Domain decomposition
% Omega_1: Square [ 1,3]x[-1,1], Cartesian coords
% Omega_2: Square [-1,1]x[-1,1], Cartesian coords
% Omega_3: Square [-3,1]x[-1,1], Cartesian coords

m([1,2,3])=m; m1=m(1); m2=m(2); m3=m(3);
E1=eye(m1);
E2=eye(m2);
E3=eye(m3);
E4=eye(n);

% Omega_1
[Dx1,x1]=chebD(m1);
x1=x1+2;
A1=Dx1*Dx1;
g1=m1; rd1=[1,m1]; kd1=2:m1-1;
a1=[1,1]; b1=[0,0];

% Omega_2
[Dx2,x2]=chebD(m2);
x2=x2+0;
A2=Dx2*Dx2;
g2=[1,m2]; rd2=[1,m2]; kd2=2:m2-1;
a2=[1,1]; b2=[0,0];

% Omega_3
[Dx3,x3]=chebD(m3);
x3=x3-2;
A3=Dx3*Dx3;
g3=1; rd3=[1,m3]; kd3=2:m3-1;
a3=[1,1]; b3=[0,0];

% Same vertical coordinate y in [-1,1]
[Dy,y]=chebD(n); y=y';
B=Dy*Dy;
rd4=[1,n]; kd4=2:n-1;
a4=[1,1];  b4=[0,0];

%% Coordinate Transformation

% Tensor product grid
[xx1,yy1]=ndgrid(x1,y);
[xx2,yy2]=ndgrid(x2,y);
[xx3,yy3]=ndgrid(x3,y);

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
F2=0*xx2;
F3=-sin(pi*xx3)./(pi*xx3);

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

%% Interior Green functions

[Lx,Ly]=ndgrid(L1,L4); LL1=1./(Lx+Ly);
[Lx,Ly]=ndgrid(L2,L4); LL2=1./(Lx+Ly);
[Lx,Ly]=ndgrid(L3,L4); LL3=1./(Lx+Ly);
GF1=@(rhs) S1*(LL1.*(S1(kd1,:)\rhs/S4(kd4,:).'))*S4.';
GF2=@(rhs) S2*(LL2.*(S2(kd2,:)\rhs/S4(kd4,:).'))*S4.';
GF3=@(rhs) S3*(LL3.*(S3(kd3,:)\rhs/S4(kd4,:).'))*S4.';
%% Schur complement method

t1  = Dx1(g1,kd1)*S1(kd1,:)*diag(S1(kd1,:)\A1(kd1,g1))*LL1;
t3  = Dx3(g3,kd3)*S3(kd3,:)*diag(S3(kd3,:)\A3(kd3,g3))*LL3;
t211=-Dx2(g2(1),kd2)*S2(kd2,:)*diag(S2(kd2,:)\A2(kd2,g2(1)))*LL2;
t222=-Dx2(g2(2),kd2)*S2(kd2,:)*diag(S2(kd2,:)\A2(kd2,g2(2)))*LL2;
t212=-Dx2(g2(1),kd2)*S2(kd2,:)*diag(S2(kd2,:)\A2(kd2,g2(2)))*LL2;
t221=-Dx2(g2(2),kd2)*S2(kd2,:)*diag(S2(kd2,:)\A2(kd2,g2(1)))*LL2;


Sig11=S4(kd4,:)*diag(Dx1(g1,g1)-Dx2(g2(1),g2(1))-t1-t211)/S4(kd4,:);
Sig22=S4(kd4,:)*diag(Dx3(g3,g3)-Dx2(g2(2),g2(2))-t3-t222)/S4(kd4,:);
Sig12=S4(kd4,:)*diag(-Dx2(g2(1),g2(2))-t212)/S4(kd4,:);
Sig21=S4(kd4,:)*diag(-Dx2(g2(2),g2(1))-t221)/S4(kd4,:);
Sig=[Sig11, Sig12; Sig21, Sig22];

f1=-Dx1(g1,kd1)*K1(GF1(K1(F1)))+Dx2(g2(1),kd2)*K2(GF2(K2(F2)));
f2=-Dx3(g3,kd3)*K3(GF3(K3(F3)))+Dx2(g2(2),kd2)*K2(GF2(K2(F2)));
f=[f1,f2];

ub=zeros(2,n);
ub(:,kd4)=reshape(Sig\f(:),[],2)';
ub(:,rd4)=ub(:,kd4)*G4';

%% Subdomain solutions

uu1=GF1(K1(F1-OP1(E1(:,g1)*ub(1,:))))+E1(:,g1)*ub(1,:);
uu2=GF2(K2(F2-OP2(E2(:,g2)*ub)))+E2(:,g2)*ub;
uu3=GF3(K3(F3-OP3(E3(:,g3)*ub(2,:))))+E3(:,g3)*ub(2,:);


%% Plot Solution
figure(1); clf; hold on;

% Plot Omega_1
surf(xx1,yy1,uu1);

% Plot Omega_2
surf(xx2,yy2,uu2);

% Plot Omega_2
surf(xx3,yy3,uu3);
hold off;

axis manual;
colormap(jet(256));
shading interp; camlight;
xlabel('x'); ylabel('y');

end