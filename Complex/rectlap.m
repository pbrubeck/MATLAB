function [] = rectlap(vertex, corners, N, k)
% Dirichlet eigenmodes of the laplacian on simple polygons mapped to the
% rectangle
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
Dxx=(2/dx)^2*(Dx*Dx);
Dyy=(2/dy)^2*(Dy*Dy);

% Jacobian determinant
[xx,yy]=ndgrid(xmin+dx*(x+1)/2, ymin+dy*(y+1)/2);
zz=xx+1i*yy;
J=abs(evaldiff(f,zz(2:end-1,2:end-1))).^2;

% Compute the k-th eigenmode
[V1, w1]=eig(Dxx(2:end-1,2:end-1),'vector'); U1=inv(V1);
[V2, w2]=eig(Dyy(2:end-1,2:end-1),'vector'); U2=inv(V2);
[W1,W2]=ndgrid(w1,w2); W=W1+W2;
OP=@(x) reshape(greenF(V1,V2,U1,U2,W,J,reshape(x, size(J))), size(x));

[V,lam]=eigs(OP, numel(J), k, 'sm');
lam=diag(lam);
phi=zeros(N);
phi(2:end-1,2:end-1)=reshape(V(:,k), size(J));
phi=phi/max(abs(phi(:)));
disp(lam);

% Conformal mapping ww=f(zz)
b1=f(zz([1 end],:)); b2=f(zz(:,[1 end]));
BC=Dxx(:,[1 end])*b1+b2*Dyy(:,[1 end])';

ww=zeros(N); ww([1 end],:)=b1; ww(:,[1 end])=b2;
ww(2:end-1,2:end-1)=greenF(V1,V2,U1,U2,W,1,-BC(2:end-1,2:end-1));
uu=real(ww); vv=imag(ww);

% Plot solution
figure(1);
surfl(uu,vv,phi,'light');
view(2); shading interp; colormap(jet(256));
zrange=max(phi(:))-min(phi(:));
xrange=max(uu(:))-min(uu(:));
yrange=max(vv(:))-min(vv(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlim([min(uu(:)) max(uu(:))]);
ylim([min(vv(:)) max(vv(:))]);
title(sprintf('\\lambda_{%d} = %.8f', k, lam(k)));
end

function X=greenF(V1, V2, U1, U2, W, J, X)
X=V1*((U1*(J.*X)*U2')./W)*V2';
end