function [lam] = meshSchurNKP(z0, quads, curv, m, k)
% Schur complement of the NKP preconditioner with quadrileteral domains.
kd=2:m-1; 
rd=[1,m];
vec=@(x) x(:);

tol=1E-15;
maxit=200;
restart=7;

% Differential operators
[Dx,x0]=legD(m); 
Dx=Dx(end:-1:1,end:-1:1); 
x0=x0(end:-1:1);
% Constraint operator
a=[1,1;1,1];
b=[0,0;0,0];
C1=zeros(2,m); C1(:,rd)=eye(2);
C2=zeros(2,m); C2(:,rd)=eye(2);
C1=diag(a(1,:))*C1+diag(b(1,:))*Dx(rd,:);
C2=diag(a(2,:))*C2+diag(b(2,:))*Dx(rd,:);
% xx fine grid for over-integration
[xx,wx]=gauleg(-1,1,m); xx=xx(end:-1:1);

[net, adj, corners, edges] = meshtopo(quads);
ndom=size(quads,1);
d=[ndom*(m-2)^2, size(adj,1)*(m-2), max(corners(:))];
dofs=sum(d);

% Function handles
stiff=cell(ndom,1); % block stiffness
mass =cell(ndom,1); % block mass
nkp  =cell(ndom,1); % block NKP
gf   =cell(ndom,1); % block NKP Green's function

% Construct Schur NKP preconditioner
S=sparse(m*size(adj,1), m*size(adj,1));

% Evaluate Jacobian determinant J and metric tensor [g11, g12; g12, g22]
% Galerkin stiffness and mass (matrix-free) operators, with their NKP
% Update NKP Schur complement and compute local Green functions
% This for loop is highly parallelizable, the bottleneck are the updates to
% the Schur complement sparse matrix.
for j=1:ndom
F=curvedquad(z0(quads(j,:)),curv(j,:));
[jac,g11,g12,g22] = diffgeom(F,xx,xx);
[stiff{j},mass{j},A1,B1,A2,B2]=lapGalerkin(Dx,x0,xx,wx,jac,g11,g12,g22);
[S,nkp{j},gf{j}]=feedSchurNKP(S,net(j,:),A1,B1,A2,B2,C1,C2);
end

% Schur LU decomposition
ix=1:size(S,2);
iy=zeros(m,size(adj,1));
e=ones(size(iy));
iy(2:end-1,:)=reshape(1:d(2),m-2,[]);
edges(edges==0)=-d(2);
iy([1,end],:)=d(2)+edges';
e([1,end],:)=1/2;
iy=reshape(iy, size(ix));
e=reshape(e, size(ix));
e(iy==0)=[];
ix(iy==0)=[];
iy(iy==0)=[];
Rschur=sparse(ix,iy,e, size(S,2), d(2)+d(3));
S=Rschur'*S*Rschur;

if(size(S,1) < 12E5)
    figure(2);
    imagesc(log(abs(S)));
    title(sprintf('NKP Schur complement \\Sigma\ncond(\\Sigma) = %.3f', condest(S)));
    colormap(gray(256)); colorbar; 
    axis square;
    drawnow;
end

[Lschur,Uschur]=lu(S);

function [m1,m2]=mapToCornerIndex(id)
    m1=rd(bitand(id-1,1)==[0,1]);
    m2=rd(bitand(id-1,2)==[0,2]);
end

function [m11,m12,m21,m22]=mapToEdgeIndex(edge)
    m11=rd(edge(1)==[1,2]);
    m12=rd(edge(1)==[3,4]);     
    m21=rd(edge(3)==[1,2]);
    m22=rd(edge(3)==[3,4]);
    if numel(m11)==0
        m11=kd;
    end
    if numel(m12)==0
        m12=kd;
    end
    if numel(m21)==0
        m21=kd;
    end
    if numel(m22)==0
        m22=kd;
    end
end

function [vv] = fullop(op,uu)
    vv=reshape(uu,m,m,[]);
    for r=1:size(vv,3)
        vv(:,:,r)=op{r}(vv(:,:,r));
    end
    vv=reshape(vv,size(uu));
end

function [u] = precond(rhs)
    RHS=reshape(rhs(1:d(1)), m-2, m-2, []);
    v=zeros(m,m,size(RHS,3));
    for r=1:size(v,3)
        v(:,:,r)=gf{r}(RHS(:,:,r),0,0,0);
    end
    p=d(1);
    s1=zeros(m-2, size(adj,1));
    for r=1:size(adj,1) % Parallelizable
        [m11,m12,m21,m22]=mapToEdgeIndex(adj(r,:));
        s1(:,r)=rhs(1+p:p+m-2) + ...
                -vec(nkp{adj(r,2)}(v(:,:,adj(r,2)),m11,m12)) + ...
                -vec(nkp{adj(r,4)}(v(:,:,adj(r,4)),m21,m22));
        p=p+m-2;
    end
    s0=rhs(p+1:end);
    [m1, m2] = ndgrid(rd,rd); 
    for r=1:size(corners,1)
        for t=1:size(corners,2)
            c=corners(r,t);
            if(c>0)
                s0(c) = s0(c)-nkp{r}(v(:,:,r),m1(t),m2(t));
            end
        end
    end
    srhs=[s1(:); s0(:)];

    % Solve for boundary nodes
    bb=Uschur\(Lschur\srhs);
    b1=reshape(bb(1:d(2)), m-2, []);
    b1=[zeros(m-2,1), b1];
    b2=[0; bb(d(2)+1:end)];
    
    % Solve for interior nodes with the given BCs
    u=zeros(size(v));
    for r=1:ndom % This is the most parallelizable loop
        b0=reshape(b2(1+corners(r,:)),[2,2]);
        u(:,:,r)=gf{r}(RHS(:,:,r), b1(:,1+net(r,1:2))', b1(:,1+net(r,3:4)), b0);
    end
    u=u(:);
end

function [u] = pick(uu)
    uu=reshape(uu,m,m,[]);
    u=zeros(dofs,1);
    u(1:d(1))=uu(kd,kd,:);
    p=d(1);
    for r=1:size(adj,1)
        [m1,m2]=mapToEdgeIndex(adj(r,:));
        u(1+p:p+m-2)=vec(uu(m1,m2,adj(r,2)));
        p=p+m-2;
    end
    for r=1:max(corners(:))
        [q,c]=find(corners==r,1,'first');
        [m1,m2]=mapToCornerIndex(c);
        u(1+p)=uu(m1,m2,q);
        p=p+1;
    end
end

function [uu] = assembly(u)
    uu=zeros(m,m,ndom);
    uu(kd,kd,:)=reshape(u(1:d(1)), m-2,m-2,[]);
    ub=u(1+d(1):d(1)+d(2));
    ub=reshape(ub, m-2, []);
    ub=[zeros(m-2,1), ub];
    uc=u(1+d(1)+d(2):end);
    uc=[0;uc];
    for r=1:size(uu,3)
        uu(rd,kd,r)=ub(:,1+net(r,1:2))';
        uu(kd,rd,r)=ub(:,1+net(r,3:4));
        uu(rd,rd,r)=reshape(uc(1+corners(r,:)),[2,2]);
    end
end

function [u] = Rtransp(uu)
    uu=reshape(uu,m,m,[]);
    u=zeros(dofs,1);
    u(1:d(1))=uu(kd,kd,:);
    p=d(1);
    for r=1:size(adj,1)
        [m11,m12,m21,m22]=mapToEdgeIndex(adj(r,:));
        u(1+p:p+m-2)=vec(uu(m11,m12,adj(r,2)))+vec(uu(m21,m22,adj(r,4)));
        p=p+m-2;
    end
    for r=1:max(corners(:))
        [q,c]=find(corners==r);
        for t=1:numel(c)
            [m1,m2]=mapToCornerIndex(c(t));
            u(1+p)=u(1+p)+uu(m1,m2,q(t));
        end
        p=p+1;
    end
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

[xx,yy]=ndgrid(x0);

if nargin>4
    tol=1E-11;
    [U,lam,~,~,relres]=lobpcg(rand(dofs,k),@afun,@bfun,@pfun,[],tol,maxit);
    relres=relres(:,end);
    display(relres);
    
    uuu=zeros(m,m,ndom,k);
    for j=1:k
        uuu(:,:,:,j)=assembly(U(:,j));
    end
    uuu=ipermute(uuu,[1,2,4,3]);
    
    figure(1); zoom off; pan off; rotate3d off;
    for j=1:size(uuu,4)
        mapping=curvedquad(z0(quads(j,:)),curv(j,:));
        ww=mapping(xx,yy);
        modegallery(real(ww), imag(ww), uuu(:,:,:,j));
        if j==1, hold on; end
    end
    hold off;    
else
    F=ones(m,m,ndom);
    [uu,~,relres,~]=poissonSolver(F);
    display(relres);
    lam=[];
    
    % Testing preconditioner
    %uu=reshape(precond(Rtransp(fullop(mass,F))),size(F));  
    
    figure(1);
    for j=1:ndom
        mapping=curvedquad(z0(quads(j,:)),curv(j,:));
        ww=mapping(xx,yy);
        surf(real(ww), imag(ww), uu(:,:,j));
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