function [ ] = quadmeshdemo( m )
% Attempting to combine Schur complement and NKP with quadrileteral
% domains.

% Differential operators
n=m;
[Dx,x]=chebD(m);
[Dy,y]=chebD(n);
% Constraint operator, Dirichlet BCs
C1=zeros(2,m); C1([1,end],[1,end])=eye(2);
C2=zeros(2,n); C2([1,end],[1,end])=eye(2);
% (xx,yy) fine grid for over-integration
[xx,wx]=gaulob(-1,1,m+32);
[yy,wy]=gaulob(-1,1,n+32);
% (xxx,yyy) is two dimensional, here we evaluate the equation coefficients
[xxx,yyy]=ndgrid(xx,yy);

% Vertices
v0=[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)];
% v0=2/(2-sqrt(2))*[1i;0;1];

% Sides
L=abs(v0([3,1,2])-v0([2,3,1]));
% Contact triangle
V=eye(3)+diag((sum(L)/2-L)./L([2,3,1]))*[-1,0,1; 1,-1,0; 0,1,-1];

z0=zeros(7,1);
z0([1,2,3])=v0;       % Vertices
z0([4,5,6])=V*v0;     % Touch points
z0(7)=(L'*v0)/sum(L); % Incenter

% Assemble quads [EN, ES; WN, WS]
Z1=reshape(z0([7,4,5,1]),[2,2]);
Z2=reshape(z0([7,5,6,2]),[2,2]);
Z3=reshape(z0([7,6,4,3]),[2,2]);

% Topology
east=1; west=2; north=3; south=4;
adj=zeros(3,4);
adj(1:3,[east,north])=[1,2; 2,3; 3,1];
T=topo(adj);

% Evaluate Jacobian determinant J and metric tensor [E, F; F, G]
[~,J1,E1,F1,G1]=mapquad(Z1,xxx,yyy);
[~,J2,E2,F2,G2]=mapquad(Z2,xxx,yyy);
[~,J3,E3,F3,G3]=mapquad(Z3,xxx,yyy);

% Construct Schur NKP preconditioner
S=sparse((m-2)*size(adj,1), (m-2)*size(adj,1));

% ( Init Schur )

% ...

% Galerkin stiffness and mass (matrix-free) operators, with their NKP
% Update NKP Schur complement and compute local Green functions

[stiff1,mass1,A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J1,E1,F1,G1);
[S,gf1,kd,gb]=feedSchurNKP(S, T(1,:), A1, B1, A2, B2, C1, C2);

[stiff2,mass2,A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J2,E2,F2,G2);
[S,gf2, ~, ~]=feedSchurNKP(S, T(2,:), A1, B1, A2, B2, C1, C2);

[stiff3,mass3,A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J3,E3,F3,G3);
[S,gf3, ~, ~]=feedSchurNKP(S, T(3,:), A1, B1, A2, B2, C1, C2);

figure(2);
subplot(2,2,1); imagesc(log(abs(A1))); colormap(gray(256));
subplot(2,2,2); imagesc(log(abs(B1))); colormap(gray(256));
subplot(2,2,3); imagesc(log(abs(A2))); colormap(gray(256));
subplot(2,2,4); imagesc(log(abs(B2))); colormap(gray(256));

% Helper routine to solve generalized Sylvester equations
b1=zeros(2,n);
b2=zeros(m,2);
F=ones(m,n);

afun=@(uu)  kd(stiff1(gb(uu)));
pfun=@(rhs) gf1(rhs);

ub=zeros(m,n);
rhs=kd(mass1(F)-stiff1(ub));
x0=gf1(rhs)+kd(ub);
tol=1e-13;
maxit=25;

% Solve the big problem in one domain
[uu,~,res,its]=gmres(afun,rhs,10,tol,maxit,pfun,[],x0);
uu=gb(uu)+ub;

display(its);
display(res);

[xx,yy]=ndgrid(x,y);
ww=mapquad(Z1,xx,yy);

figure(3);
surf(real(ww), imag(ww), reshape(uu,size(ww)));
colormap(jet(256));
shading interp; camlight; view(2);
dx=diff(xlim());
dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
end


function [T]=topo(adj)
% Inverts adjacency map from inter->doms to dom->inters (E,W,N,S)
    ndoms=max(adj(:));
    T=zeros(ndoms,4);
    T(adj(adj(:,1)>0,1),1)=find(adj(:,1)>0);
    T(adj(adj(:,2)>0,2),2)=find(adj(:,2)>0);
    T(adj(adj(:,3)>0,3),3)=find(adj(:,3)>0);
    T(adj(adj(:,4)>0,4),4)=find(adj(:,4)>0);
end


function [S, gf, kd, gb] = feedSchurNKP(S, T, A1, B1, A2, B2, C1, C2)
% Updates SchurNKP S given domain topology T (EWNS) and NKP (A1,B1,A2,B2)
m=size(A1,1);
n=size(B2,1);
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];

% Eigenfunctions
[K1,M1,E1,V1,L1]=eigenfunctions(A1, A2, C1);
[K2,M2,E2,V2,L2]=eigenfunctions(B2, B1, C2);

% Green's function
[Lx,Ly]=ndgrid(L1,L2);
LL=1./(Lx+Ly);
LL((Lx==0)|(Ly==0))=0;
function uu=greenF(F)
    uu=reshape(V1(kd1,kd1)*(LL.*(V1(kd1,kd1).'*reshape(F,[m-2,n-2])*V2(kd2,kd2)))*V2(kd2,kd2).',size(F));
end
gf=@greenF;
kd=@(uu) reshape(E1(:,kd1)'*uu*E2(:,kd2),[],1);
gb=@(uu) E1(:,kd1)*reshape(uu,m-2,n-2)*E2(:,kd2)';

% Topology
[I1,J1]=meshgrid(rd1(T(1:2)>0));
[I2,J2]=meshgrid(rd2(T(3:4)>0));

% Building blocks go here


% ...


% Assembly

% Schur complement block-indices
[ii,jj]=meshgrid(T(T>0));

end


function [SK,SM,E,V,L]=eigenfunctions(K,M,C)
% Computes stiffness SK and mass SM matrices in the constrained basis E and
% constrainded eigenfunctions V and eigenvalues L
m=size(K,1);
kd=2:m-1;
rd=[1,m];

% Constrained basis
E=eye(m);
E(rd,kd)=-C(:,kd);
E(rd,:)=C(:,rd)\E(rd,:);
SK=E'*K*E;
SM=E'*M*E;

% Eigenfunctions
V=zeros(m);
[V(kd,kd),L]=eig(SK(kd,kd), SM(kd,kd), 'vector');
%V(kd,kd)=bsxfun(@rdivide, V(kd,kd), sqrt(diag(V(kd,kd)'*SM(kd,kd)*V(kd,kd)))');
V(rd,kd)=E(rd,kd)*V(kd,kd);
end