function [q] = HelmConformal(a, b, N, k)
f=sqrt(a^2-b^2);
xi0=acosh(a/f);

[Dx,x]=chebD(2*N-1); x=xi0*x; Dx=Dx/xi0; Dxx=Dx*Dx;
[A1,G]=setBC(Dxx,Dx,1,0);
y=2*pi/N*(1:N);
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=f*cosh(zz);
uu=real(ww);
vv=imag(ww);
J=abs(f*sinh(zz)).^2;
J=J(2:end-1,:);

% Factor homogeneous BC Green function
[V1, L1]=eig(A1,'vector'); W1=inv(V1);
L2=-[0:N/2-1 -N/2:-1].^2;
[L1,L2]=ndgrid(L1,L2); LL=L1+L2;
function X=Afun(X)
    X=reshape(X, size(J));
    X=real(V1*fft((W1*ifft(J.*X,[],2))./LL,[],2));
    X=X(:);
end

% Compute k eigenmodes
[V,lam]=eigs(@Afun, numel(J), k, 'sm');
[lam, id]=sort(-diag(lam));
q=lam.^2*f^2/4;
psi=zeros(size(ww));
psi(2:end-1,:)=reshape(V(:,id(k)), size(J));
psi([1,end],:)=G*psi(2:end-1,:);
psi=psi/max(abs(psi(:)));


% Plot solution
figure(1);
surf(uu(2:N,[1:end,1]),vv(2:N,[1:end,1]),psi(2:N,[1:end,1]));
camlight; view(2); colormap(jet(256)); shading interp;

zrange=max(psi(:))-min(psi(:));
xrange=max(uu(:))-min(uu(:));
yrange=max(vv(:))-min(vv(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlim([min(uu(:)) max(uu(:))]);
ylim([min(vv(:)) max(vv(:))]);
title(sprintf('q_{%d} = %.8f', k, q(k)));
end