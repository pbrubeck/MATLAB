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

% Jabian determinant
[xx,yy]=ndgrid(xmin+dx*(x+1)/2, ymin+dy*(y+1)/2);
zz=xx+1i*yy;
J=abs(evaldiff(f,zz(2:end-1,2:end-1))).^2;

% Compute the k-th eigenmode
[V1, w1]=eig(Dxx,'vector'); U1=inv(V1);
[V2, w2]=eig(Dyy,'vector'); U2=inv(V2);
[W1,W2]=ndgrid(w1,w2); 
W=W1+W2;
OP=@(x) reshape(invop(V1,V2,U1,U2,W,J,reshape(x, size(J))), size(x));
[V,lam]=eigs(OP, numel(J), k, 'sm');
lam=diag(lam);
[lam,idx]=sort(lam,'descend');
phi=zeros(N);
phi(2:end-1,2:end-1)=reshape(V(:,idx(k)), N-2);
disp(lam);

% Conformal mapping
ww=f(zz);
uu=real(ww);
vv=imag(ww);

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
title(sprintf('\\lambda_{%d} = %.8f', k, lam(k)));
end

function X=invop(V1, V2, U1, U2, W, J, X)
X=V1*((U1*(J.*X)*V2)./W)*U2;
end