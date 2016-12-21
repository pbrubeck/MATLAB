function [lam] = rectlap(vertex, corners, N, k)
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
[Dx,x]=chebD(N(1)); Dx=2/dx*Dx;
[Dy,y]=chebD(N(2)); Dy=2/dy*Dy;
A1=Dx*Dx;
A2=Dy*Dy;

% Jacobian determinant
[xx,yy]=ndgrid(xmin+dx*(x+1)/2, ymin+dy*(y+1)/2);
zz=xx+1i*yy;
J=abs(evaldiff(f,zz(2:end-1,2:end-1))).^2;

% Compute the k-th eigenmode
[V1, L1]=eig(A1(2:end-1,2:end-1),'vector'); U1=inv(V1);
[V2, L2]=eig(A2(2:end-1,2:end-1),'vector'); U2=inv(V2);
[L1,L2]=ndgrid(L1,L2); LL=L1+L2;
OP=@(x) reshape(poissonSquare(V1,V2,U1,U2,LL,J,reshape(x, size(J))), size(x));

[V,lam]=eigs(OP, numel(J), k, 'sm');
lam=diag(lam);
phi=zeros(N);
phi(2:end-1,2:end-1)=reshape(V(:,k), size(J));
phi=phi/max(abs(phi(:)));

% Conformal mapping ww=f(zz)
b1=f(zz([1 end],:)); b2=f(zz(:,[1 end]));
BC=A1(:,[1 end])*b1+b2*A2(:,[1 end])';

ww=zeros(N); ww([1 end],:)=b1; ww(:,[1 end])=b2;
ww(2:end-1,2:end-1)=poissonSquare(V1,V2,U1,U2,LL,1,-BC(2:end-1,2:end-1));
uu=real(ww); vv=imag(ww);

% Plot solution
surf(uu,vv,phi);
colormap(jet(256)); view(2); shading interp;
title(sprintf('\\lambda_{%d}=%.8f', k, lam(k)));
xlim([min(uu(:)) max(uu(:))]);
ylim([min(vv(:)) max(vv(:))]);
xlabel('x'); ylabel('y'); axis square; axis manual;
end

function X=poissonSquare(V1, V2, U1, U2, LL, J, X)
X=V1*((U1*(J.*X)*U2')./LL)*V2';
end