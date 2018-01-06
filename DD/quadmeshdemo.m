function [ ] = quadmeshdemo( m )
% Attempting to combine Schur complement and NKP with quadrileteral
% domains.
n=m;
kd1=2:m-1;
kd2=2:n-1;
vec=@(x) x(:);

% Differential operators
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
v0=2/(2-sqrt(2))*[1i;0;1]; %Right angle
v0=[-2+3i;0;5]; % Scalene
v0=[3i;-1;1]; % Isoceles
v0=[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)]; % Equilateral

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
adj(1:3,[east,north])=[2,1; 3,2; 1,3];
net=topo(adj);

% Evaluate Jacobian determinant J and metric tensor [E, F; F, G]
[~,J1,E1,F1,G1]=mapquad(Z1,xxx,yyy);
[~,J2,E2,F2,G2]=mapquad(Z2,xxx,yyy);
[~,J3,E3,F3,G3]=mapquad(Z3,xxx,yyy);

% Construct Schur NKP preconditioner
S=sparse((m-2)*size(adj,1), (m-2)*size(adj,1));

% Galerkin stiffness and mass (matrix-free) operators, with their NKP
% Update NKP Schur complement and compute local Green functions
stiff=cell(3,1);
mass =cell(3,1);
[stiff{1},mass{1},A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J1,E1,F1,G1);
[S,nkp1,gf1]=feedSchurNKP(S, net(1,:), A1, B1, A2, B2, C1, C2);

[stiff{2},mass{2},A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J2,E2,F2,G2);
[S,nkp2,gf2]=feedSchurNKP(S, net(2,:), A1, B1, A2, B2, C1, C2);

[stiff{3},mass{3},A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J3,E3,F3,G3);
[S,nkp3,gf3]=feedSchurNKP(S, net(3,:), A1, B1, A2, B2, C1, C2);

% Full operators
function [v]=fullmass(u)
    v=reshape(u,m,n,[]);
    for i=1:size(v,3)
        v(:,:,i)=mass{i}(v(:,:,i));
    end
    v=reshape(v,size(u));
end
function [v]=fullstiff(u)
    v=reshape(u,m,n,[]);
    for i=1:size(v,3)
        v(:,:,i)=stiff{i}(v(:,:,i));
    end
    v=reshape(v,size(u));
end

% Schur LU decomposition
[Lschur, Uschur, pschur]=lu(S,'vector');

figure(2);
imagesc(log(abs(S)));
title(sprintf('cond(\\Sigma) = %f', condest(S)));
colormap(gray(256)); colorbar; axis square;
drawnow;

net(net==0)=max(net(:))+1;

% Schur NKP Preconditioner (Triangle exclusive, working on generalization)
function [u] = precond(F)
    F=reshape(F, m, n, []);
    v=zeros(m,n,3);
    v(:,:,1)=gf1(F(:,:,1),0,0);
    v(:,:,2)=gf2(F(:,:,2),0,0);
    v(:,:,3)=gf3(F(:,:,3),0,0);
    adj
    net
    srhs=zeros(m-2, size(adj,1));
    srhs(:,1)=vec(F(1,kd2,adj(1,1))-nkp1(v(:,:,adj(1,1)),1,kd2)) + ...
              vec(F(kd1,1,adj(1,3))-nkp2(v(:,:,adj(1,3)),kd1,1));
    srhs(:,2)=vec(F(1,kd2,adj(2,1))-nkp2(v(:,:,adj(2,1)),1,kd2)) + ...
              vec(F(kd1,1,adj(2,3))-nkp3(v(:,:,adj(2,3)),kd1,1));
    srhs(:,3)=vec(F(1,kd2,adj(3,1))-nkp3(v(:,:,adj(3,1)),1,kd2)) + ...
              vec(F(kd1,1,adj(3,3))-nkp1(v(:,:,adj(3,3)),kd1,1));
    srhs=srhs(:);
    
    % Solve for boundary nodes
    b=Uschur\(Lschur\srhs(pschur));
    b=reshape(b, m-2, []);
    b=[b, zeros(m-2,1)];
    
    % Solve for interior nodes with the given BCs
    % Not ready, should add extra params to gf
    u=zeros(size(F));
    u(:,:,1)=gf1(F(:,:,1), b(:,net(1,1:2))', b(:,net(1,3:4)));
    u(:,:,2)=gf2(F(:,:,2), b(:,net(2,1:2))', b(:,net(2,3:4)));
    u(:,:,3)=gf3(F(:,:,3), b(:,net(3,1:2))', b(:,net(3,3:4)));
    
    u=u(:);
end

% Right-hand sides
F=ones(m,n,3);
ub=zeros(m,n,3);
rhs=fullmass(F)-fullstiff(ub);
uu=reshape(precond(rhs), size(rhs));


% % Helper routine to solve generalized Sylvester equations
% afun=@(uu)  kd(stiff{1}(gb(uu)));
% pfun=@(rhs) kd(gf1(gb(rhs),0,0));
% 
% F=ones(m,n);
% ub=zeros(m,n);
% rhs=kd(mass{1}(F)-stiff{1}(ub));
% x0=kd(gf1(mass{1}(F)-stiff{1}(ub),0,0))+kd(ub);
% tol=1e-14;
% maxit=25;
% 
% % Solve the big problem in one domain
% [uu,~,res,its]=gmres(afun,rhs,5,tol,maxit,pfun,[],x0);
% uu=gb(uu)+ub;
% 
% display(its);
% display(res);


[xx,yy]=ndgrid(x,y);
ww1=mapquad(Z1,xx,yy);
ww2=mapquad(Z2,xx,yy);
ww3=mapquad(Z3,xx,yy);

figure(1);
surf(real(ww1), imag(ww1), uu(:,:,1)); hold on;
surf(real(ww2), imag(ww2), uu(:,:,2));
surf(real(ww3), imag(ww3), uu(:,:,3)); hold off;

colormap(jet(256));
shading interp; camlight; view(2);
dx=diff(xlim());
dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
end

function [net]=topo(adj)
% Inverts adjacency map from inter->doms to dom->inters (E,W,N,S)
ndoms=max(adj(:));
net=zeros(ndoms,4);
net(adj(adj(:,1)>0,1),1)=find(adj(:,1)>0);
net(adj(adj(:,2)>0,2),2)=find(adj(:,2)>0);
net(adj(adj(:,3)>0,3),3)=find(adj(:,3)>0);
net(adj(adj(:,4)>0,4),4)=find(adj(:,4)>0);
end