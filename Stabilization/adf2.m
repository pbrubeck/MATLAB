function [relres,it,resvec,icolor] = adf2(n,ne,nu,opts)
if(nargin==3)
    opts=ones(1,5);
end
opts=[opts(:)',ones(1,5-numel(opts))];
no=opts(1);
ifras=logical(opts(2));
ifcrs=logical(opts(3));
ifneu=logical(opts(4));
ifdeal=logical(opts(5));
ifsweep=no<0;
no=min(max(0,no),1);

ifplot=true;

% GMRES settings
maxit=200;
tol=1e-10;

% Problem settings
if(nu<0), nu=1./abs(nu); end

%UDATA=@(x,y) x.*(1-exp((y-1)/nu))./(1-exp(-2/nu)); FDATA=0; bcbox=[1,1,1,1];
%UDATA=0; FDATA=1; bcbox=[1,1,1,0];
UDATA=1; FDATA=0; bcbox=[1,1,1,1];
ux=@(x,y) 0+2*y.*(1-x.^2);
uy=@(x,y) 0-2*x.*(1-y.^2);

% Stabilization CFL
CFL=1.0E-2;
%CFL=inf;

function [u]=vel(x,y,ijac)
    u=ijac(:,:,:,:,1).*ux(x,y)+ijac(:,:,:,:,2).*uy(x,y);
end


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

x=zhat*hx+repmat(conv(xc,[1,1]/2,'valid'),n,1);
y=zhat*hy+repmat(conv(yc,[1,1]/2,'valid'),n,1);
[x,y]=ndgrid(x,y);
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
vx=ux(xm(:),ym(:));
vy=uy(xm(:),ym(:));

%--------------------------------------------------------------------------
%  Boundary conditions
%--------------------------------------------------------------------------

bc=zeros(nfaces,nel);
bc(:)=2; % Overlap
icolor=1:nel;

if(ifsweep)
    wx=hx';
    wy=hy';
    
    itopo=box_topo(nex,ney);
    iflux=get_graph(itopo,vx,vy,wx,wy);
        
    %ex=1; ey=1;
    ex=floor((nex+1)/2)+1; ey=ney;
    e=ex+nex*(ey-1);
    iflux(:,e)=max(0,iflux(:,e));
    [isweep,icolor]=toposort_loops(itopo,iflux);
  
    bc=reshape(bc,[],nex,ney); 
    icolor=reshape(icolor,nex,ney);
    bc(1,2:end  ,:)=1+(icolor(2:end  ,:)<=icolor(1:end-1,:)); 
    bc(2,1:end-1,:)=1+(icolor(1:end-1,:)<=icolor(2:end  ,:)); 
    bc(3,:,2:end  )=1+(icolor(:,2:end  )<=icolor(:,1:end-1)); 
    bc(4,:,1:end-1)=1+(icolor(:,1:end-1)<=icolor(:,2:end  ));    
        
    bc=reshape(bc,[],nel);
    icolor=reshape(icolor,1,[]);
end

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
if(ifdeal)
    nxd=ceil(3*n/2);
else
    nxd=n;
end
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

kx=repmat(1:length(xc),2,1); kx=kx(2:end-1);
ky=repmat(1:length(yc),2,1); ky=ky(2:end-1);
[xxc,yyc]=ndgrid(xc(kx),yc(ky),1);
quad=reshape(permute(reshape(xxc+1i*yyc,[2,nex,2,ney]),[1,3,2,4]),[4,nex*ney]);

function [map]=mymap(e,x,y,z)
    rad=[inf;inf;inf;inf];
    w=crvquad(quad([4,3,2,1],e),rad);
    map=[real(w(x,y)); imag(w(x,y)); z];
end

xd=zeros(nxd,nyd,nzd,nel);
yd=zeros(nxd,nyd,nzd,nel);
zd=zeros(nxd,nyd,nzd,nel);
bm1=zeros(nxd,nyd,nzd,nel);
G=zeros(nxd,nyd,nzd,6,nel);
C=zeros(nxd,nyd,nzd,3,nel);

iffast=(numel(uniquetol(vx,tol))==1)&(numel(uniquetol(vy,tol))==1)&...
       (numel(uniquetol(hx,tol))==1)&(numel(uniquetol(hy,tol))==1);

ncoef=nel;
if(iffast), ncoef=1; end
for e=1:ncoef
% Equation coefficients
coord=@(x,y,z) mymap(e,x,y,z);
[xd(:,:,:,e),yd(:,:,:,e),zd(:,:,:,e),ijac,bm1(:,:,:,e),G(:,:,:,:,e)]=gmfact(ndim,coord,xq,wq);
ud=vel(xd(:,:,:,e),yd(:,:,:,e),ijac);
C(:,:,:,1,e)=bm1(:,:,:,e).*ud(:,:,:,1);
C(:,:,:,2,e)=bm1(:,:,:,e).*ud(:,:,:,2);
C(:,:,:,3,e)=bm1(:,:,:,e).*ud(:,:,:,3);
end
if(iffast)
    bm1=repmat(bm1(:,:,:,1),1,1,1,nel);
    G=repmat(G(:,:,:,:,1),1,1,1,1,nel);
    C=repmat(C(:,:,:,:,1),1,1,1,1,nel);
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

function [bx]=bfun(x,elems)
    if(nargin==1)
        elems=1:nel;
    end
    elems=mod(elems-1,nel)+1;
    bx=reshape(x,n,n,nel);
    for ie=elems
        bx(:,:,ie)=J'*(bm1(:,:,:,ie).*(J*bx(:,:,ie)*J'))*J;
    end
    bx=dssum(bx);
    bx=reshape(bx,size(x));
end

function [au]=afun(u,elems)
    if(nargin==1)
        elems=1:nel;
    end
    sz=size(u);
    u=reshape(u,n,n,nel);
    au=zeros(size(u));
    elems=mod(elems-1,nel)+1;
    for ie=elems
        fu=F*u(:,:,ie)*F';
        u_x=D*u(:,:,ie)*J';
        u_y=J*u(:,:,ie)*D';
        au(:,:,ie) = D'*(G(:,:,:,1,ie).*u_x + G(:,:,:,4,ie).*u_y)*J + ...
                     J'*(G(:,:,:,4,ie).*u_x + G(:,:,:,2,ie).*u_y)*D + ...
                     J'*(C(:,:,:,1,ie).*u_x + C(:,:,:,2,ie).*u_y + ...
                     (1/dt)*bm1(:,:,:,ie).*(J*(u(:,:,ie)-fu)*J'))*J;
    end
    au=dssum(au);
    au(:)=mask(:).*au(:);
    au=reshape(au,sz);
end


%--------------------------------------------------------------------------
%  Preconditioner
%--------------------------------------------------------------------------
[A,B]=schwarz2d(n,no,hx,hy,nu,vx,vy,dt,bc,ifdeal,ifneu);

% Schwarz weight
if(ifsweep)
    wt=ones(n,n,nel);
    wt(1,:,bc(1,:)==1)=0;
    wt(n,:,bc(2,:)==1)=0;
    wt(:,1,bc(3,:)==1)=0;
    wt(:,n,bc(4,:)==1)=0;
    wt=reshape(wt,n,n,nel,[]);
elseif(ifras || no==0)
    wt=ones(n,n,nel);
elseif(no==1)
    wt1=ones(ns,ns,nel);
    wt2=ones(ns,ns,nel);
    wt1=extrude(wt1,0,0.0,wt2,0,1.0);
    wt2=dssum(wt2);
    wt2=extrude(wt2,0,1.0,wt1,0,-1.0);
    wt2=extrude(wt2,2,1.0,wt2,0,1.0);
    wt=wt2(2:end-1,2:end-1,:);
end
wt=wt./dssum(wt); % ensure partition of unity
wt(mask==0)=0;

function [v1]=extrude(v1,l1,f1,v2,l2,f2)
    k1=[1+l1,size(v1,1)-l1];
    k2=[1+l2,size(v2,1)-l2];
    v1(k1,2:end-1,:)=f1*v1(k1,2:end-1,:)+f2*v2(k2,2:end-1,:);
    k1=[1+l1,size(v1,2)-l1];
    k2=[1+l2,size(v2,2)-l2];
    v1(2:end-1,k1,:)=f1*v1(2:end-1,k1,:)+f2*v2(2:end-1,k2,:);
end

function [u]=psweep(r)
    visit=zeros(nex,ney);
    im=zeros(length(xc),length(yc));
    figure(2); hv=pcolor(xc,yc,im');
    colormap(gray(2));
    caxis('manual'); caxis([0,1]);
    set(gca,'YDir','normal');
    title('Sweeping');  
    
    u=zeros(size(r));
    w=zeros(size(r));
    ncolor=max(icolor);
    for ic=1:ncolor
        ie=find(icolor==ic);
        w(:,:,ie)=r(:,:,ie)-w(:,:,ie);
        z=pschwarz(w,ie);
        u(:,:,ie)=u(:,:,ie)+z(:,:,ie);
        u=wt.*u;
        u=dssum(u);
        if(ic<ncolor)
        je=find(icolor==ic|icolor==ic+1|icolor==ic+2);
        w=afun(u,je);
        end    
        visit(ie)=visit(ie)+1;
        im(1:nex,1:ney)=visit;
        set(hv,'CData',im');
        drawnow;
    end
end

function [u]=pschwarz(r,elems)
    if(nargin==1)
        elems=1:size(r,3);
    end
    u=reshape(r,size(A,1),size(B,1),[]);
    for ie=elems
        je=mod(ie-1,size(u,3))+1;
        LA=A(:,:,2,ie)\A(:,:,1,ie);
        LB=B(:,:,2,ie)\B(:,:,1,ie);
        LR=A(:,:,2,ie)\u(:,:,je)/B(:,:,2,ie)';
        u(:,:,je)=sylvester(LA,LB',LR);
%         K=kron(B(:,:,2,ie),A(:,:,1,ie))+kron(B(:,:,1,ie),A(:,:,2,ie))+...
%           kron(B(:,:,3,ie),A(:,:,3,ie));        
%         u(:,:,je)=reshape(K\reshape(u(:,:,je),[],1),size(u(:,:,ie)));
    end
    u=reshape(u,size(r));
end

function [u]=psmooth(r)
    if(ifsweep)
    u=psweep(reshape(r,n,n,[]));
    elseif(no==0)
    u=pschwarz(reshape(r,n,n,[]));
    elseif(no==1)
    % go to exteded array
    w1=zeros(ns,ns,nel);
    w1(2:end-1,2:end-1,:)=reshape(r,n,n,[]);
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
    end
    % sum border nodes    
    u(:)=wt(:).*u(:);
    u=dssum(u);
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


function [u]=pfun(r)
    u=psmooth(r);
    if(ifcrs)
        u=u+pcoarse(r-afun(u));
    end
end

%--------------------------------------------------------------------------

mult=1./dssum(ones(n,n,nel));

ub=zeros(n,n,nel);
f=zeros(n,n,nel);
if(isfloat(UDATA))
    ub(x==1)=UDATA;
else
    ub(mask==0)=UDATA(x(mask==0),y(mask==0));
end
if(isfloat(FDATA))
    f(:)=FDATA;
else
    f=FDATA(x(mask==0),y(mask==0));
end
b=bfun(f)-afun(ub);

u0=pcoarse(b(:));
u0(:)=0;

restart=maxit;
[u,flag,relres,it,resvec]=gmres(@afun,b(:),restart,tol,1,[],@pfun,u0(:));
% relres=0; resvec=0; %u=pfun(b);
it=length(resvec)-1;
if(flag>0)
    it=maxit;
end

for k=1:0
u0=u0+pfun(b(:)-afun(u0));
end
%u=u0;


u=reshape(u,size(ub));
u=u+ub;
%u=UDATA(x,y)-reshape(u,size(x));


figure(2);
semilogy(0:length(resvec)-1,resvec*relres/resvec(end),'.-b');
ylim([tol/100,1]);
drawnow;

if(ifplot)
figure(1);
semplot2(x,y,u); 
%shading interp; camlight;
%view(2);
colormap(jet(256));
colorbar;
end

return;
if(n*n*nel<=(5*8)^2)
    figure(7);
    gfun = @(x) afun(pfun(x))-x;
    fullcond(gfun,n,nex,ney);
    hold on;
    plot(exp(1i*linspace(0,2*pi,512)),'k');
    hold off;
    axis equal;
    title(sprintf('Spectrum of P^{-1}A-I, Pe=%1.1E',1/nu));
end
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
    surf(x(:,:,e),y(:,:,e),u(:,:,e)); hold on; %drawnow;
end
hold off; 
end

function [kap]=fullcond(afun,n,nex,ney)
    
    gidx = repmat((1:n)',1,nex)+(n-1)*repmat(0:nex-1,n,1);
    gidy = repmat((1:n)',1,ney)+(n-1)*repmat(0:ney-1,n,1);
    
    
    [gx,gy]=ndgrid(gidx(:),gidy(:));
    gid=gx+((n-1)*nex+1)*(gy-1);
    gid=reshape(permute(reshape(gid,n,nex,n,ney),[1,3,2,4]),n,n,nex,ney);
    
    mask=ones(size(gid));
    mask(1,:,1,:)=0;
    %mask(n,:,nex,:)=0;
    mask(:,1,:,1)=0;
    %mask(:,n,:,ney)=0;

    ntot=max(gid(:));
    x1=zeros(ntot,1);
    y1=zeros(ntot,1);

    amat=zeros(ntot,ntot);
    for j=1:ntot
        x1(j)=1;
        y1(gid)=afun(x1(gid));
        amat(:,j)=y1(:);
        x1(j)=0;
    end

    x1(gid)=mask;
    p=x1==1;
    amat=amat(p,p);
    
    lam=eig(amat);
    plot(real(lam),imag(lam),'.b');
    kap=cond(amat);
    %kap=0;
end