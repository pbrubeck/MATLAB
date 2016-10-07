function [] = DirichletEllipse(a,b,N,M)
% Solves Dirichlet BVP over an elliptical region
f=sqrt(a^2-b^2);
xi0=acosh(a/f);

[Dx,x]=chebD(2*N-1); x=xi0*x; Dx=Dx/xi0; Dxx=Dx*Dx;
[Dyy,y]=fourD2(M);
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=f*cosh(zz);
uu=real(ww);
vv=imag(ww);
J=abs(f*sinh(zz)).^2;

% Boundary conditions
bc=[cos(5*y);cos(5*y)];
F=zeros(size(zz));
RHS=J.*F-Dxx(:,[1,end])*bc;

psi=zeros(size(zz));
psi([1 end],:)=bc;
psi(2:end-1,:)=sylvester(Dxx(2:end-1,2:end-1), Dyy', RHS(2:end-1,:));

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
end