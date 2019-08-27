function [relres,it,resvec] = adf2(n,ne,nu)

bcbox=[1,0,1,0];

ux=1;
uy=1;
CFL=0.02;
%CFL=inf;

maxit=200;
tol=1e-10;

no=1;
ns=n+2*no;
ndim=2;
nfaces=2*ndim;

nex=ne(1);
ney=ne(end);
nel=nex*ney;

% SEM hat
[Dhat,zhat,what]=legD(n);
Bhat=diag(what);
Ahat=Dhat'*Bhat*Dhat;
Chat=Bhat*Dhat;
Jfem=[1-zhat, 1+zhat]/2;
F=bubfilt(zhat);

xc=linspace(-1,1,nex+1);
yc=linspace(-1,1,ney+1);
hx=diff(xc)/2;
hy=diff(yc)/2;

x1=zhat*hx+repmat(conv(xc,[1,1]/2,'valid'),n,1);
y1=zhat*hy+repmat(conv(yc,[1,1]/2,'valid'),n,1);
[x,y]=ndgrid(x1,y1);
x=semreshape(x,n,nex,ney);
y=semreshape(y,n,nex,ney);

dx=min(diff(zhat))*sqrt(min(kron(hx,hy)));
dt=CFL*dx/max(hypot(ux,uy),[],'all');

[hx,hy]=ndgrid(hx,hy);
hx=hx(:);
hy=hy(:);

vx=zeros(nel,1);
vy=zeros(nel,1);
vx(:)=ux;
vy(:)=uy;

%--------------------------------------------------------------------------
%  Boundary conditions
%--------------------------------------------------------------------------

bc=zeros(nfaces,nel);
bc(:)=2; % Overlap

% Override with problem boundaries
bc=reshape(bc,nfaces,nex,ney);
bc(1,  1,:)=bcbox(1); % Left
bc(2,end,:)=bcbox(2); % Right
bc(3,:,1  )=bcbox(3); % Bottom
bc(4,:,end)=bcbox(4); % Top
bc=reshape(bc,nfaces,nel);

mask=ones(n*nex,n*ney);
mask(1  ,:)=mask(1  ,:)*(1-bcbox(1));
mask(end,:)=mask(end,:)*(1-bcbox(2));
mask(:,1  )=mask(:,1  )*(1-bcbox(3));
mask(:,end)=mask(:,end)*(1-bcbox(4));
mask=semreshape(mask,n,nex,ney);


%--------------------------------------------------------------------------
%  Forward operator
%--------------------------------------------------------------------------

nxd=ceil(3*n/2); %nxd=n;
nyd=nxd; nzd=1;
if(nxd==n)
    [xq,wq]=deal(zhat,what);
    J=eye(n);
else
    [xq,wq]=gauleg(-1,1,nxd);
    J=legC(zhat,xq);
end
D=J*Dhat;
wq=reshape(kron(wq,wq),[nxd,nxd]);

kx=repmat(1:length(xc),2,1);
ky=repmat(1:length(yc),2,1);
kx=kx(2:end-1);
ky=ky(2:end-1);
[xxc,yyc]=ndgrid(xc(kx),yc(ky),1);
vert=reshape(permute(reshape(xxc+1i*yyc,[2,nex,2,ney]),[1,3,2,4]),[4,nex*ney]);

function [map]=mymap(e,x,y,z)
    rad=[inf;inf;inf;inf];
    w=crvquad(vert([4,3,2,1],e),rad);
    map=[real(w(x,y)); imag(w(x,y)); z];
end

xd=zeros(nxd,nyd,nzd,nel);
yd=zeros(nxd,nyd,nzd,nel);
zd=zeros(nxd,nyd,nzd,nel);
bm1=zeros(nxd,nyd,nzd,nel);
G=zeros(nxd,nyd,nzd,6,nel);
C=zeros(nxd,nyd,nzd,3,nel);
for e=1:nel
% Equation coefficients
coord=@(x,y,z) mymap(e,x,y,z);
[xd(:,:,:,e),yd(:,:,:,e),zd(:,:,:,e),ijac,bm1(:,:,:,e),G(:,:,:,:,e)]=gmfact(ndim,coord,xq,wq);
C(:,:,:,1,e)=bm1(:,:,:,e).*ux;
C(:,:,:,2,e)=bm1(:,:,:,e).*uy;
end
G=nu.*G;

function [u]=dssum(u)
    sz=size(u);
    u=reshape(u,size(u,1),size(u,2),nex,ney);
    s=u(1,:,2:end,:)+u(end,:,1:end-1,:);
    u(1,:,2:end,:)=s;
    u(end,:,1:end-1,:)=s;
    s=u(:,1,:,2:end)+u(:,end,:,1:end-1);
    u(:,1,:,2:end)=s;
    u(:,end,:,1:end-1)=s;
    u=reshape(u,sz);
end

function [bx]=bfun(x)
    pe=1;
    for q=1:ney
        for p=1:nex
            x(:,:,p,q)=J'*(bm1(:,:,:,pe).*(J*x(:,:,p,q)*J'))*J;
            pe=pe+1;
        end
    end
    bx=dssum(x);
end

function [au]=afun(u,ids)
    if(nargin==1)
        ids=1;
    end
    sz=size(u);
    u=reshape(u,n,n,nex,ney);
    au=zeros(size(u));
    pe=1;
    for q=1:ney
        for p=1:nex
        fu=F*u(:,:,p,q)*F';
        ux=D*u(:,:,p,q)*J';
        uy=J*u(:,:,p,q)*D';
        au(:,:,p,q) = D'*(G(:,:,:,1,pe).*ux + G(:,:,:,4,pe).*uy)*J + ...
                      J'*(G(:,:,:,4,pe).*ux + G(:,:,:,2,pe).*uy)*D + ...
                      J'*(C(:,:,:,1,pe).*ux + C(:,:,:,2,pe).*uy + ...
                      (1/dt)*bm1(:,:,:,pe).*(J*(u(:,:,p,q)-fu)*J'))*J;
        pe = pe+1;
        end
    end
    if(ids==1)
        au=dssum(au);
        au(:)=mask(:).*au(:);
    end
    au=reshape(au,sz);
end


%--------------------------------------------------------------------------
%  Preconditioner
%--------------------------------------------------------------------------

A=zeros(ns,ns,2,nel);
B=zeros(ns,ns,2,nel);
for e=1:nel
[A(:,:,:,e),B(:,:,:,e)]=schwarz2d(n,no,hx(e),hy(e),nu,vx(e),vy(e),CFL,bc(:,e));
end

wt1=ones(ns,ns,nel);
wt2=ones(ns,ns,nel);
wt1=extrude(wt1,0,0.0,wt2,0,1.0);
wt2=dssum(wt2);
wt2=extrude(wt2,0,1.0,wt1,0,-1.0);
wt2=extrude(wt2,2,1.0,wt2,0,1.0);
wt=wt2(2:end-1,2:end-1,:);
wt=dssum(wt);
wt(:)=mask(:)./wt(:);

function [v1]=extrude(v1,l1,f1,v2,l2,f2)
    k1=[1+l1,size(v1,1)-l1];
    k2=[1+l2,size(v2,1)-l2];
    v1(k1,2:end-1,:)=f1*v1(k1,2:end-1,:)+f2*v2(k2,2:end-1,:);
    k1=[1+l1,size(v1,2)-l1];
    k2=[1+l2,size(v2,2)-l2];
    v1(2:end-1,k1,:)=f1*v1(2:end-1,k1,:)+f2*v2(2:end-1,k2,:);
end

function [u]=pschwarz(r)
    u=zeros(size(r));
    for ie=1:size(u,3)
        LA=A(:,:,2,ie)\A(:,:,1,ie);
        LB=B(:,:,2,ie)\B(:,:,1,ie);
        LR=A(:,:,2,ie)\r(:,:,ie)/B(:,:,2,ie)';
        u(:,:,ie)=sylvester(LA,LB',LR);
    end
end

function [u]=psmooth(r)
    % go to exteded array
    w1=zeros(ns,ns,nel);
    w1(2:end-1,2:end-1,:)=reshape(wt(:).*r(:),n,n,[]);
    % exchange interior nodes
    w1=extrude(w1,0,0.0,w1,2,1.0);
    w1=dssum(w1);
    w1=extrude(w1,0,1.0,w1,2,-1.0);
    % do the local solves
    w2=pschwarz(w1);
    % sum overlap region
    w1=extrude(w1,0,0.0,w2,0,1.0);
    w2=dssum(w2);
    w2=extrude(w2,0,1.0,w1,0,-1.0);
    w2=extrude(w2,2,1.0,w2,0,1.0);
    % go back to regular size array
    u=w2(2:end-1,2:end-1,:);
    % sum border nodes
    u=dssum(u);
    u=reshape(u,size(r));
end

function [u]=pcoarse(r)
    u=zeros(size(r));
end

function [u]=pfun(r)
    u=psmooth(r);
end

%--------------------------------------------------------------------------


f=1*ones(n,n,nex,ney);
ub=zeros(n,n,nex,ney);
b=bfun(f);
b=bfun(f)-afun(ub);


restart=maxit;
[u,flag,relres,it,resvec] = gmres(@afun,b(:),restart,tol,1,@pfun,[]);
%relres=0; resvec=0; u=pfun(b);
it=length(resvec)-1;

u=reshape(u,size(ub));
u=u+ub;
figure(1);
semplot2(x,y,u); 
%shading interp; camlight;
%view(2);
colormap(jet(256));
colorbar;
end


function [x]=semreshape(x,n,nex,ney)
    x=permute(reshape(x,n,nex,n,ney),[1,3,2,4]);
end

function []=semplot2(x,y,u)
u=reshape(u,size(u,1),size(u,2),[]);
x=reshape(x,size(u));
y=reshape(y,size(u));
nel = size(u,3);
for e=1:nel
    surf(x(:,:,e),y(:,:,e),u(:,:,e));hold on;
end
hold off; 
end