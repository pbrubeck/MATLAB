function [lam] = testSchurNKP( m, k )
% Schur complement of the NKP preconditioner with quadrileteral domains.
n=m;
kd1=2:m-1; rd1=[1,m];
kd2=2:n-1; rd2=[1,n];
vec=@(x) x(:);

tol=1E-15;
maxit=100;
restart=7;

% Differential operators
[Dx,x0]=chebD(m);
[Dy,y0]=chebD(n);
% Constraint operator
a=[1,1;1,1];
b=[0,0;0,0];
C1=zeros(2,m); C1(:,rd1)=eye(2);
C2=zeros(2,n); C2(:,rd2)=eye(2);
C1=diag(a(1,:))*C1+diag(b(1,:))*Dx(rd1,:);
C2=diag(a(2,:))*C2+diag(b(2,:))*Dy(rd2,:);
% (xx,yy) fine grid for over-integration
[xx,wx]=gaulob(-1,1,m+32);
[yy,wy]=gaulob(-1,1,n+32);
% (xxx,yyy) is two dimensional, here we evaluate the equation coefficients
[xxx,yyy]=ndgrid(xx,yy);



% Vertices
v0=[-2+3i;0;2];                                           % Scalene
v0=4/sqrt(3*sqrt(3))*[1i;exp(1i*pi*7/6);exp(-1i*pi*1/6)]; % Equilateral
v0=[2i;-1;1];                                             % Isoceles
v0=2/(2-sqrt(2))*[1i;0;1];                                % Right angle



L=abs(v0([3,1,2])-v0([2,3,1])); % Sides
V=eye(3)+diag((sum(L)/2-L)./L([2,3,1]))*[-1,0,1; 1,-1,0; 0,1,-1];

z0=zeros(7,1);
z0([1,2,3])=v0;       % Vertices
z0([4,5,6])=V*v0;     % Touch points
z0(7)=(L'*v0)/sum(L); % Incenter

% Assemble quads [EN, ES; WN, WS]
Z=zeros(2,2,3);
Z(:,:,1)=z0([7,4;5,1]);
Z(:,:,2)=z0([7,5;6,2]);
Z(:,:,3)=z0([7,6;4,3]);

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
for j=1:ndom
[~,jac,g11,g12,g22]=mapquad(Z(:,:,j),xxx,yyy);
[stiff{j},mass{j},A1,B1,A2,B2]=lapGalerkin(Dx,Dy,x0,y0,xx,yy,wx,wy,jac,g11,g12,g22);
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
    RHS=reshape(rhs(1:d(1)), m-2, n-2, []);
    v=zeros(m,n,size(RHS,3));
    for r=1:size(adj,1)
        v(:,:,r)=gf{r}(RHS(:,:,r),0,0,0);
    end
    p=d(1);
    s1=zeros(m-2, size(adj,1));
    for r=1:size(adj,1)
        s1(:,r)=rhs(1+p:p+m-2) + ...
                -vec(nkp{adj(r,1)}(v(:,:,adj(r,1)),1,kd2)) + ...
                -vec(nkp{adj(r,3)}(v(:,:,adj(r,3)),kd1,1));
        p=p+m-2;
    end
    s0=rhs(p+1)-nkp{1}(v(:,:,1),1,1)-nkp{2}(v(:,:,2),1,1)-nkp{3}(v(:,:,3),1,1);
    srhs=[s1(:); s0];

    % Solve for boundary nodes
    bb=Uschur\(Lschur\srhs(pschur));
    b1=reshape(bb(1:end-corners), m-2, []);
    b1=[b1, zeros(m-2,1)];
    b0=zeros(2,2); 
    b0(1,1)=bb(end-corners+1:end);

    % Solve for interior nodes with the given BCs
    u=zeros(size(v));
    for r=1:ndom
        u(:,:,r)=gf{r}(RHS(:,:,r), b1(:,net(r,1:2))', b1(:,net(r,3:4)), b0);
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
        u(1+p:p+m-2)=vec(uu(bx,kd2,adj(r,1)));
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
        u(1+p:p+m-2)=vec(uu(bx,kd2,adj(r,1)))+vec(uu(kd1,by,adj(r,3)));
        p=p+m-2;
    end
    u(1+p)=sum(uu(1,1,:));
end

function [u] = afun(u)
    for r=1:size(u,2)
        u(:,r)=Rtransp(fullop(stiff, assembly(u(:,r))));
    end
end

function [u] = bfun(u)
    for r=1:size(u,2)
        u(:,r)=Rtransp(fullop(mass, assembly(u(:,r))));
    end
end

pcalls=0;
function [u] = pfun(u)
    for r=1:size(u,2)
        pcalls=pcalls+1;
        u(:,r)=pick(precond(u(:,r)));
    end    
end

function [uu,flag,relres,iter]=poissonSolver(F,ub)
    if nargin==1
        ub=zeros(size(F));
    end
    rhs=Rtransp(fullop(mass,F)-fullop(stiff,ub));
    uu=pfun(rhs);
    [uu,flag,relres,iter]=gmres(@afun,rhs,restart,tol,ceil(maxit/restart),[],@pfun,uu);
    uu=ub+reshape(assembly(uu),size(ub));  
end

Nq=256;
xq=linspace(-1,1,Nq); %xq=x;
yq=linspace(-1,1,Nq); %yq=y;
[xx,yy]=ndgrid(xq,yq);

if nargin>1
    tol=1E-11;
    [U,lam,~,~,relres]=lobpcg(rand(dofs,k),@afun,@bfun,@pfun,[],tol,maxit);
    uuu=zeros(m,n,ndom,k);
    for j=1:k
        uuu(:,:,:,j)=assembly(U(:,j));
    end
    uuu=ipermute(uuu,[1,2,4,3]);
    
    figure(1); zoom off; pan off; rotate3d off;
    for j=1:size(uuu,4)
        ww=mapquad(Z(:,:,j),xx,yy);
        modegallery(real(ww), imag(ww), interpcheb(interpcheb(uuu(:,:,:,j),xq,1),yq,2));
        if j==1, hold on; end
    end
    hold off;    
else
    F=ones(m,n,ndom);
    [uu,~,relres,~]=poissonSolver(F);
    lam=[];
    display(relres);
    
    figure(1);
    for j=1:ndom
        ww=mapquad(Z(:,:,j),xx,yy);
        surf(real(ww), imag(ww), interpcheb(interpcheb(uu(:,:,j),xq,1),yq,2));
        if j==1, hold on; end
    end 
    hold off;
end

colormap(jet(256));
view(2);
shading interp; camlight;
dx=diff(xlim());
dy=diff(ylim());
pbaspect([dx,dy,min(dx,dy)]);
display(pcalls);
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