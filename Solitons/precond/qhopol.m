function [] = qhopol(m,n)
L=10;
omega=2/10;
lam=0;
VL=@(r) (omega*r).^2; 

[rr,th,jac,M,H,U]=schrodpol(m,n,L,lam,VL);
xx=rr.*cos(th);
yy=rr.*sin(th);
ii=1:m;
jj=[1:n,1];

nx=0;
ny=0;
u0 = HermitePsi([zeros(nx,1);1],xx*sqrt(omega)).*...
     HermitePsi([zeros(ny,1);1],yy*sqrt(omega));
u=u0/sqrt(M(u0,u0));

t=0;
setlatex();
mytitle='$z = %f$, $E/\\omega = %f$, $P = %f$';
figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L,L]);
ylim([-L,L]);
axis square;
shading interp;
colormap(magma(256));
colorbar();
view(2);
E=real(H(u,u))/(2*omega);
P=real(M(u,u));
title(sprintf(mytitle,t,E,P));
drawnow;

nframes=0;
T=2*pi/omega;
dt=T/nframes;
for k=1:nframes
    u=U(dt,u);
    t=t+dt;
    set(h1,'ZData',abs(u(ii,jj)).^2);
    E=real(H(u,u))/(2*omega);
    P=real(M(u,u));
    title(sprintf(mytitle,t,E,P));
    drawnow;
end
end