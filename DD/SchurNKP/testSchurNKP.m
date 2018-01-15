function [ ] = testSchurNKP( m )
% Schur complement of the NKP preconditioner with quadrileteral domains.
n=m;
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];
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
v0=[3i;-1;1];                                             % Isoceles
v0=[-2+3i;0;2];                                           % Scalene
v0=2/(2-sqrt(2))*[1i;0;1];                                % Right angle
v0=2/sqrt(3*sqrt(3))*[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)]; % Equilateral


L=abs(v0([3,1,2])-v0([2,3,1])); % Sides
V=eye(3)+diag((sum(L)/2-L)./L([2,3,1]))*[-1,0,1; 1,-1,0; 0,1,-1];

z0=zeros(7,1);
z0([1,2,3])=v0;       % Vertices
z0([4,5,6])=V*v0;     % Touch points
z0(7)=(L'*v0)/sum(L); % Incenter

% Assemble quads [EN, ES; WN, WS]
Z=zeros(2,2,3);
Z(:,:,1)=reshape(z0([7,4,5,1]),[2,2]);
Z(:,:,2)=reshape(z0([7,5,6,2]),[2,2]);
Z(:,:,3)=reshape(z0([7,6,4,3]),[2,2]);

% Topology
adj=zeros(3,4);
adj(1:3,[1,3])=[2,1; 3,2; 1,3]; % [E,W,N,S]
net=topo(adj);
corners=zeros(size(net));
corners(:,1)=1;                 % [EN,WN,ES,WS]
net=[net, corners];
ndom=max(adj(:));

% Evaluate Jacobian determinant J and metric tensor [E, F; F, G]
J=zeros([size(xxx),ndom]);
E=zeros([size(xxx),ndom]);
F=zeros([size(xxx),ndom]);
G=zeros([size(xxx),ndom]);
f=cell(3,1);
for j=1:ndom
    [~,J(:,:,j),E(:,:,j),F(:,:,j),G(:,:,j)]=mapquad(Z(:,:,j),xxx,yyy);
%     [f{j},J(:,:,j)]=confmapquad(Z(:,:,j),xxx,yyy);
%     J(:,:,j)=-J(:,:,j);
%     E(:,:,j)=J(:,:,j);
%     G(:,:,j)=J(:,:,j);
end

% Construct Schur NKP preconditioner
S11=sparse((m-2)*size(adj,1),(m-2)*size(adj,1));
S12=sparse((m-2)*size(adj,1), max(corners(:)) );
S21=sparse( max(corners(:)) ,(m-2)*size(adj,1));
S22=sparse( max(corners(:)) , max(corners(:)) );

% Galerkin stiffness and mass (matrix-free) operators, with their NKP
% Update NKP Schur complement and compute local Green functions
stiff=cell(3,1);
mass =cell(3,1);
nkp  =cell(3,1);
gf   =cell(3,1);
for j=1:3
[stiff{j},mass{j},A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,J(:,:,j),E(:,:,j),F(:,:,j),G(:,:,j));
[S11,S12,S21,S22,nkp{j},gf{j}]=feedSchurNKP(S11,S12,S21,S22,net(j,:),A1,B1,A2,B2,C1,C2);
end

% Schur LU decomposition
S=[S11,S12;S21,S22];
[Lschur, Uschur, pschur]=lu(S,'vector');

figure(2);
imagesc(log(abs(S)));
title(sprintf('cond(\\Sigma) = %.3f\ncond(\\Sigma_{11}) = %.3f', condest(S), condest(S11)));
colormap(gray(256)); colorbar; axis square;
drawnow;

net(net==0)=max(net(:))+1; 

function [vv]=fullop(op,uu)
    vv=reshape(uu,m,n,[]);
    for r=1:size(vv,3)
        vv(:,:,r)=op{r}(vv(:,:,r));
    end
    vv=reshape(vv,size(uu));
end

function [u] = precond(F)
    F=reshape(F, m, n, []);
    v=zeros(size(F));
    for r=1:size(adj,1)
        v(:,:,r)=gf{r}(F(:,:,r),0,0,0);
    end
    
    s1=zeros(m-2, size(adj,1));
    for r=1:size(adj,1)
        s1(:,r)=vec(F(1,kd2,adj(r,1))-nkp{adj(r,1)}(v(:,:,adj(r,1)),1,kd2))+...
                vec(F(kd1,1,adj(r,3))-nkp{adj(r,3)}(v(:,:,adj(r,3)),kd1,1));
    end
    s0=F(1,1,1)-nkp{1}(v(:,:,1),1,1)+...
       F(1,1,2)-nkp{2}(v(:,:,2),1,1)+...
       F(1,1,3)-nkp{3}(v(:,:,3),1,1);
    srhs=[s1(:); s0];

    % Solve for boundary nodes
    b=Uschur\(Lschur\srhs(pschur));
    b1=reshape(b(1:end-corners), m-2, []);
    b1=[b1, zeros(m-2,1)];
    b0=zeros(2,2); 
    b0(1,1)=b(end-corners+1:end);

    % Solve for interior nodes with the given BCs
    u=zeros(size(F));
    for r=1:size(adj,1)
        u(:,:,r)=gf{r}(F(:,:,r), b1(:,net(r,1:2))', b1(:,net(r,3:4)), b0);
    end
    u=u(:);
end

function [b]=afun(uu)
    g=max(net(:));
    uu=reshape(uu, m,n,[]);
    bc=zeros(size(uu));
    for r=1:size(uu,3)
        bx=rd1(net(r,1:2)==g);
        by=rd2(net(r,3:4)==g);
        bc(bx,kd2,r)=uu(bx,kd2,r);
        bc(kd1,by,r)=uu(kd1,by,r);
        bc(bx,by,r)=uu(bx,by,r);
        uu(bx,kd2,r)=0;
        uu(kd1,by,r)=0;
        uu(bx,by,r)=0;
    end
    vv=reshape(fullop(stiff, uu), m,n,[]);
    for r=1:size(uu,3)
        bx=rd1(net(r,1:2)==g);
        by=rd2(net(r,3:4)==g);
        vv(bx,kd2,r)=bc(bx,kd2,r);
        vv(kd1,by,r)=bc(kd1,by,r);
        vv(bx,by,r)=bc(bx,by,r);   
    end
    b=vv(:);
end

function [b]=pfun(rhs)
    g=max(net(:));
    rhs=reshape(rhs, m,n,[]);
    vv=reshape(precond(rhs), m,n,[]);
    for r=1:size(rhs,3)
        bx=rd1(net(r,1:2)==g);
        by=rd2(net(r,3:4)==g);
        vv(bx,kd2,r)=rhs(bx,kd2,r);
        vv(kd1,by,r)=rhs(kd1,by,r);
        vv(bx,by,r)=rhs(bx,by,r);    
    end
    b=vv(:);
end

% Right-hand sides
F=ones(m,n,ndom);
ub=zeros(m,n,ndom);
rhs=fullop(mass,F)-fullop(stiff,ub);
uu=precond(rhs);

tol=1e-14;
maxit=25;

% [uu,~,res,its]=gmres(@afun,rhs(:),5,tol,maxit,@pfun,[],uu);
% display(its);
% display(res);
uu=ub+reshape(uu,m,n,[]);

[xx,yy]=ndgrid(x,y);
figure(1); clf; grid on; hold on;
for j=1:ndom
    ww=mapquad(Z(:,:,j),xx,yy);
%     ww=gf{j}(zeros(m,n), ...
%         f{j}(xx(rd1,kd2),yy(rd1,kd2)), ...
%         f{j}(xx(kd1,rd2),yy(kd1,rd2)), ...
%         f{j}(xx(rd1,rd2),yy(rd1,rd2)));
    surf(real(ww), imag(ww), uu(:,:,j));
end
hold off;

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