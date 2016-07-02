function [] = rectlap(vertex, corners, N, k)
% Dirichlet eigenmodes of the laplacian on simple polygons
p=polygon(vertex);
f=rectmap(p, corners);
params=parameters(f);
xmin=min(real(params.prevertex));
xmax=max(real(params.prevertex));
ymin=min(imag(params.prevertex));
ymax=max(imag(params.prevertex));
dx=xmax-xmin;
dy=ymax-ymin;

% Set differential operators
N([1,2])=N;
[Dx,x]=chebD(N(1));
[Dy,y]=chebD(N(2));
Dxx=(2/dx)^2*Dx*Dx;
Dxx=Dxx(2:end-1,2:end-1);
Dyy=(2/dy)^2*Dy*Dy;
Dyy=Dyy(2:end-1,2:end-1);

% Conformal mapping
[xx,yy]=ndgrid(xmin+dx*(x+1)/2, ymin+dy*(y+1)/2);
zz=xx+1i*yy;
ww=f(zz);
uu=real(ww);
vv=imag(ww);
df=evaldiff(f,zz);
J=abs(df).^2;

phi=zeros(N);
%B=-ones(N);
%RHS=J.*B;
%phi(2:end-1,2:end-1)=sylvester(Dxx,Dyy',RHS(2:end-1,2:end-1));

% Compute the k-th eigenmode
N2=(N(1)-2)*(N(2)-2);
L=kron(speye(size(Dyy,1)),sparse(Dxx))+kron(sparse(Dyy),speye(size(Dxx,1)));

J=J(2:end-1,2:end-1);
J=sparse(1:N2,1:N2,J(:));
[V,lam]=eigs(L, J, k,'sm');
lam=diag(lam);
[lam,idx]=sort(lam,'descend');
V=V(:,idx);
phi(2:end-1,2:end-1)=reshape(V(:,k), [N(1)-2, N(2)-2]);

% Plot solution
figure(1);
surfl(uu,vv,phi,'light'); 
shading interp;
colormap(jet(256));
zrange=max(phi(:))-min(phi(:));
xrange=max(uu(:))-min(uu(:));
yrange=max(vv(:))-min(vv(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlim([min(uu(:)) max(uu(:))]);
ylim([min(vv(:)) max(vv(:))]);
view(0,90);
title(sprintf('\\lambda_{%d} = %f', k, lam(k)));
end