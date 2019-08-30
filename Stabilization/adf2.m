function [relres,it,resvec] = adf2(n,ne,nu)
% GMRES settings
maxit=100;
tol=1e-10;

% Problem settings
bcbox=[1,0,1,1];
UDATA=0;
FDATA=1;

% Velocity field
ux=@(x,y,z) 1+0*y.*(1-x.^2);
uy=@(x,y,z) 0-0*x.*(1-y.^2);
uz=@(x,y,z) 0+0*x+0*y;

% Stabilization CFL
CFL=1.0E-2;
%CFL=inf;

function [u]=vel(x,y,z,ijac)
    u=ijac(:,:,:,:,1).*ux(x,y,z)+ijac(:,:,:,:,2).*uy(x,y,z)+ijac(:,:,:,:,3).*uz(x,y,z);
end
if(nu<0)
    nu=1./abs(nu);
end

ifras=true;
no=1;
ns=n+2*no;
ndim=2;
nfaces=2*ndim;

nex=ne(1);
ney=ne(end);
nel=nex*ney;

[Dhat,zhat,what]=legD(n);
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

dx=min(diff(zhat))*min([hx(:);hy(:)]);
dt=CFL*dx;

[hx,hy]=ndgrid(hx,hy);
hx=hx(:);
hy=hy(:);

xm=conv(xc,[1,1]/2,'valid');
ym=conv(yc,[1,1]/2,'valid');
[xm,ym]=ndgrid(xm,ym);
zm=zeros(size(xm));
vx=ux(xm(:),ym(:),zm(:));
vy=uy(xm(:),ym(:),zm(:));


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

% Mask
mask=ones(n*nex,n*ney);
mask(1  ,:)=mask(1  ,:)*(1-bcbox(1));
mask(end,:)=mask(end,:)*(1-bcbox(2));
mask(:,1  )=mask(:,1  )*(1-bcbox(3));
mask(:,end)=mask(:,end)*(1-bcbox(4));
mask=semreshape(mask,n,nex,ney);


%--------------------------------------------------------------------------
%  Coarse grid
%--------------------------------------------------------------------------

%bcbox=[1,0,0,0]; % Override BCs

cmask=ones(2,2,nex,ney);
cbnd=false(nex+1,ney+1);
if(bcbox(1)==1)
    cmask(1,:,1,:)=0;
    cbnd(1,:)=true;
end
if(bcbox(2)==1)
    cmask(end,:,end,:)=0;
    cbnd(end,:)=true;
end
if(bcbox(3)==1)
    cmask(:,1,:,1)=0;
    cbnd(:,1)=true;
end
if(bcbox(4)==1)
    cmask(:,end,:,end)=0;
    cbnd(:,end)=true;
end
cbnd=cbnd(:);

gidx=repmat((1:2)',1,nex)+repmat(0:nex-1,2,1);
gidy=repmat((1:2)',1,ney)+repmat(0:ney-1,2,1);
[gx,gy]=ndgrid(gidx(:),gidy(:));
gid=gx+((2-1)*nex+1)*(gy-1);
cid=reshape(permute(reshape(gid,2,nex,2,ney),[1,3,2,4]),2,2,nel);

[Acrs,supg] = crs2(cid,hx,hy,nu,vx,vy);
ibnd=(size(Acrs,1)+1)*(find(cbnd)-1)+1;
Acrs(cbnd,:)=0;
Acrs(:,cbnd)=0;
Acrs(ibnd)=1;

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
pts=reshape(permute(reshape(xxc+1i*yyc,[2,nex,2,ney]),[1,3,2,4]),[4,nex*ney]);

function [map]=mymap(e,x,y,z)
    rad=[inf;inf;inf;inf];
    w=crvquad(pts([4,3,2,1],e),rad);
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
ud=vel(xd(:,:,:,e),yd(:,:,:,e),zd(:,:,:,e),ijac);
C(:,:,:,1,e)=bm1(:,:,:,e).*ud(:,:,:,1);
C(:,:,:,2,e)=bm1(:,:,:,e).*ud(:,:,:,2);
C(:,:,:,3,e)=bm1(:,:,:,e).*ud(:,:,:,3);
end
G=nu.*G;

mult=1./dssum(ones(n,n,nel));

function [d]=semdot(a,b)
    d=(a(:).*mult(:))'*b(:);
end

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
    bx=reshape(x,n,n,nel);
    for ie=1:nel
        bx(:,:,ie)=J'*(bm1(:,:,:,ie).*(J*bx(:,:,ie)*J'))*J;
    end
    bx=dssum(bx);
    bx=reshape(bx,size(x));
end

function [au]=afun(u,ids)
    if(nargin==1)
        ids=1;
    end
    sz=size(u);
    u=reshape(u,n,n,nel);
    au=zeros(size(u));
    for ie=1:nel
        fu=F*u(:,:,ie)*F';
        u_x=D*u(:,:,ie)*J';
        u_y=J*u(:,:,ie)*D';
        au(:,:,ie) = D'*(G(:,:,:,1,ie).*u_x + G(:,:,:,4,ie).*u_y)*J + ...
                     J'*(G(:,:,:,4,ie).*u_x + G(:,:,:,2,ie).*u_y)*D + ...
                     J'*(C(:,:,:,1,ie).*u_x + C(:,:,:,2,ie).*u_y + ...
                     (1/dt)*bm1(:,:,:,ie).*(J*(u(:,:,ie)-fu)*J'))*J;
    end
    au=dssum(au);
    if(ids==1)
        au(:)=mask(:).*au(:);
    end
    au=reshape(au,sz);
end

function [fu]=ffun(u)
    fu=reshape(u,n,n,nel);
    for ie=1:nel
        fu(:,:,ie)=F*fu(:,:,ie)*F';
    end
    fu=reshape(fu,size(u));
end

%--------------------------------------------------------------------------
%  Preconditioner
%--------------------------------------------------------------------------
[A,B]=schwarz2d(n,no,hx,hy,nu,vx,vy,dt,bc);

% Schwarz weight
if(ifras)
    wt=mask;
elseif(no==1)
    wt1=ones(ns,ns,nel);
    wt2=ones(ns,ns,nel);
    wt1=extrude(wt1,0,0.0,wt2,0,1.0);
    wt2=dssum(wt2);
    wt2=extrude(wt2,0,1.0,wt1,0,-1.0);
    wt2=extrude(wt2,2,1.0,wt2,0,1.0);
    wt=wt2(2:end-1,2:end-1,:);
    wt=dssum(wt);
    wt(:)=mask(:)./wt(:);
elseif(no==0)
    wt=reshape(mask(:).*mult(:),n,n,nel);
end

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
    if(no==1)
    % go to exteded array
    w1=zeros(ns,ns,nel);
    w1(2:end-1,2:end-1,:)=reshape(wt(:).*r(:),n,n,[]);
    % exchange interior nodes
    w1=extrude(w1,0,0.0,w1,2,1.0);
    w1=dssum(w1);
    w1=extrude(w1,0,1.0,w1,2,-1.0);
    % do the local solves
    w2=pschwarz(w1); 
    if(~ifras)
    % sum overlap region
    w1=extrude(w1,0,0.0,w2,0,1.0);
    w2=dssum(w2);
    w2=extrude(w2,0,1.0,w1,0,-1.0);
    w2=extrude(w2,2,1.0,w2,0,1.0);
    end
    % go back to regular size array
    u=w2(2:end-1,2:end-1,:);
    % sum border nodes
    u=dssum(u);
    elseif(no==0)
    u=dssum(pschwarz(wt.*reshape(r,size(wt))));
    end
    if(ifras)
    u(:)=mult(:).*u(:);
    end
    u=reshape(u,size(r));
end

function [u]=pcoarse(r)
    sz=size(r);
    r=reshape(r(:).*mult(:),n,n,[]);
    rc=zeros(2,2,size(r,3));
    for ie=1:size(r,3)
        rc(:,:,ie)=Jfem'*r(:,:,ie)*Jfem;
    end
    rc(:)=cmask(:).*rc(:);
    uc=supg(Acrs,rc);
    uc(:)=cmask(:).*uc(:);
    u=zeros(size(r));
    for ie=1:size(r,3)
        u(:,:,ie)=Jfem*uc(:,:,ie)*Jfem';
    end
    u=reshape(u,sz);
end

sigma=0;
function [u]=pfun(r)
    u=psmooth(r);
    u=u+pcoarse(r-afun(u));
    u(:)=mask(:).*u(:);
    %return;

    sigma=2/3;
    w=r-afun(u);
    z=psmooth(w);
    if(sigma==0)
        az=afun(z);
        sigma=semdot(az,w)/semdot(az,az);
        sigma=min(1,sigma);
        display(sigma);
    end
    u=u+sigma*z;
end

%--------------------------------------------------------------------------


f=FDATA*ones(n,n,nex,ney);
ub=zeros(n,n,nex,ney);
ub(x==1)=UDATA;
b=bfun(f)-afun(ub);

u0=pcoarse(b(:));
u0(:)=0;


% restart=10;
% [u1,flag,relres,it,resvec]=gmres(@afun,b(:),restart,tol,1,[],@pfun,u0(:));


restart=maxit;
[u,flag,relres,it,resvec]=gmres(@afun,b(:),restart,tol,1,[],@pfun,u0(:));
% relres=0; resvec=0; u=pfun(b);
it=length(resvec)-1;

%u=u-pfun(b(:));
u=reshape(u,size(ub));
u=u+ub;


figure(1);
semplot2(x,y,u); 
%shading interp; camlight;
%view(2);
colormap(jet(256));
colorbar;

figure(2);
semilogy(0:it,resvec*relres/resvec(end),'.-b');
ylim([tol/100,1]);
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