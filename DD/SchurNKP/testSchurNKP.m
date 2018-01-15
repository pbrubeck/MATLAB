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
v0=[-2+3i;0;2];                                           % Scalene
v0=4/sqrt(3*sqrt(3))*[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)]; % Equilateral
v0=2/(2-sqrt(2))*[1i;0;1];                                % Right angle
v0=[2i;-1;1];                                             % Isoceles

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
d=[ndom*(m-2)^2, size(adj,1)*(m-2), max(corners(:))];
dofs=sum(d);

% Function handles
stiff=cell(ndom,1); % block stiffness
mass =cell(ndom,1); % block mass
nkp  =cell(ndom,1); % block NKP
gf   =cell(ndom,1); % block NKP Green's function

% Construct Schur NKP preconditioner
S11=sparse(d(2),d(2));
S12=sparse(d(2),d(3));
S21=sparse(d(3),d(2));
S22=sparse(d(3),d(3));

% Evaluate Jacobian determinant J and metric tensor [E, F; F, G]
% Galerkin stiffness and mass (matrix-free) operators, with their NKP
% Update NKP Schur complement and compute local Green functions
for j=1:3
[~,jac,g11,g12,g22]=mapquad(Z(:,:,j),xxx,yyy);
[stiff{j},mass{j},A1,B1,A2,B2]=lapGalerkin(Dx,Dy,xx,yy,wx,wy,jac,g11,g12,g22);
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

function [vv] = fullop(op,uu)
    vv=reshape(uu,m,n,[]);
    for r=1:size(vv,3)
        vv(:,:,r)=op{r}(vv(:,:,r));
    end
    vv=reshape(vv,size(uu));
end

function [u] = precond(rhs)
    F=reshape(rhs(1:d(1)), m-2, n-2, []);
    v=zeros(m,n,size(F,3));
    for r=1:size(adj,1)
        v(:,:,r)=gf{r}(F(:,:,r),0,0,0);
    end
    p=d(1);
    s1=zeros(m-2, size(adj,1));
    for r=1:size(adj,1)
        s1(:,r)=rhs(1+p:p+m-2) + ...
                -vec(nkp{adj(r,1)}(v(:,:,adj(r,1)),1,kd2))+ ...
                -vec(nkp{adj(r,3)}(v(:,:,adj(r,3)),kd1,1));
        p=p+m-2;
    end
    s0=rhs(p+1)-nkp{1}(v(:,:,1),1,1)-nkp{2}(v(:,:,2),1,1)-nkp{3}(v(:,:,3),1,1);
    srhs=[s1(:); s0];

    % Solve for boundary nodes
    b=Uschur\(Lschur\srhs(pschur));
    b1=reshape(b(1:end-corners), m-2, []);
    b1=[b1, zeros(m-2,1)];
    b0=zeros(2,2); 
    b0(1,1)=b(end-corners+1:end);

    % Solve for interior nodes with the given BCs
    u=zeros(size(v));
    for r=1:size(adj,1)
        u(:,:,r)=gf{r}(F(:,:,r), b1(:,net(r,1:2))', b1(:,net(r,3:4)), b0);
    end
    u=u(:);
end

function [u] = pick(uu)
    uu=reshape(uu,m,n,[]);
    u=zeros(dofs,1);
    u(1:d(1))=uu(kd1,kd2,:);
    p=d(1);
    for r=1:size(adj,1)
        bx=rd1(adj(r,1:2)>0);
        u(1+p:p+m-2)=uu(bx,kd2,adj(r,1))';
        p=p+m-2;
    end
    u(1+p)=uu(1,1,1);
end

function [uu] = assembly(u)
    uu=zeros(m,n,ndom);
    p=ndom*(m-2)^2;
    uu(kd1,kd2,:)=reshape(u(1:p), m-2,n-2,[]);
    for r=1:size(adj,1)
        bx=rd1(adj(r,1:2)>0);
        by=rd2(adj(r,3:4)>0);
        uu(bx,kd2,adj(r,1))=reshape(u(1+p:p+m-2),n-2,[])';
        uu(kd1,by,adj(r,3))=reshape(u(1+p:p+m-2),m-2,[]);
        p=p+m-2;
    end
    uu(1,1,:)=u(1+p);
end

function [u] = Rtransp(uu)
    uu=reshape(uu,m,n,[]);
    u=zeros(dofs,1);
    u(1:d(1))=uu(kd1,kd2,:);
    p=d(1);
    for r=1:size(adj,1)
        bx=rd1(adj(r,1:2)>0);
        by=rd2(adj(r,3:4)>0);
        u(1+p:p+m-2)=uu(bx,kd2,adj(r,1))'+uu(kd1,by,adj(r,3));
        p=p+m-2;
    end
    u(1+p)=sum(uu(1,1,:));
end

function [b] = afun(u)
    b=Rtransp(fullop(stiff, assembly(u)));
end
pcalls=0;
function [u] = pfun(b)
    pcalls=pcalls+1;
    u=pick(precond(b));
end

% Right-hand sides
F=ones(m,n,ndom);
ub=zeros(m,n,ndom);
rhs=Rtransp(fullop(mass,F)-fullop(stiff,ub));
uu=pfun(rhs);

tol=1e-14;
maxit=10;
[uu,~,res,its]=gmres(@afun,rhs,7,tol,maxit,[],@pfun,uu);
display(pcalls);
display(res);
uu=ub+assembly(uu);

[xx,yy]=ndgrid(x,y);
figure(1);
for j=1:ndom
    ww=mapquad(Z(:,:,j),xx,yy);
    surf(real(ww), imag(ww), uu(:,:,j));hold on;
end 
hold off;

colormap(jet(256));
shading interp; camlight;
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