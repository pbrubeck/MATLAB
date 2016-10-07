function [q] = HelmConformal(a, b, N, k)
f=sqrt(a^2-b^2);
xi0=acosh(a/f);

[Dx,x]=chebD(2*N-1); x=xi0*x; Dx=Dx/xi0; Dxx=Dx*Dx;
[Dyy,y]=fourD2(N);
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=f*cosh(zz);
uu=real(ww);
vv=imag(ww);
J=abs(f*sinh(zz)).^2;
J=J(2:end-1,:);

% Factor homogeneous BC Green function
[V1, w1]=eig(Dxx(2:end-1,2:end-1),'vector'); U1=inv(V1);
[V2, w2]=eig(Dyy,'vector'); U2=inv(V2);
[W1,W2]=ndgrid(w1,w2); W=W1+W2;
OP=@(x) reshape(greenF(V1,V2,U1,U2,W,J,reshape(x, size(J))), size(x));

% Compute k eigenmodes
[V,lam]=eigs(OP, numel(J), k, 'sm');
[lam, id]=sort(-diag(lam));
q=lam.^2*f^2/4;
psi=zeros(size(ww));
psi(2:end-1,:)=reshape(V(:,id(k)), size(J));
psi=psi/max(abs(psi(:)));

% Plot solution
figure(1);
surfl(uu(1:N,[1:end,1]),vv(1:N,[1:end,1]),psi(1:N,[1:end,1]),'light');

view(2); shading interp; colormap(jet(256));
zrange=max(psi(:))-min(psi(:));
xrange=max(uu(:))-min(uu(:));
yrange=max(vv(:))-min(vv(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlim([min(uu(:)) max(uu(:))]);
ylim([min(vv(:)) max(vv(:))]);
title(sprintf('q_{%d} = %.8f', k, q(k)));
end

function X=greenF(V1, V2, U1, U2, W, J, X)
X=V1*((U1*(J.*X)*U2')./W)*V2';
end