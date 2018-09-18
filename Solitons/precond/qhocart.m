function [] = qhocart(m)
L=10;
omega=2/L;
lam=0;
VX=@(x) (omega*x).^2; 
VY=@(y) (omega*y).^2; 

[xx,yy,jac,M,H,U]=schrodcart(m,L,lam,VX,VY);

nx=1;
ny=2;
u0 = HermitePsi([zeros(nx,1);1],xx*sqrt(omega)).*...
     HermitePsi([zeros(ny,1);1],yy*sqrt(omega));
u=u0/sqrt(M(u0,u0));

t=0;
setlatex();
mytitle='$z = %f$, $E/\\omega = %f$, $P = %f$';
figure(1);
h1=surf(xx,yy,abs(u).^2);
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

nframes=100;
T=2*pi/omega;
dt=T/nframes;
for k=1:nframes
    u=U(dt,u);
    t=t+dt;
    set(h1,'ZData',abs(u).^2);
    E=real(H(u,u))/(2*omega);
    P=real(M(u,u));
    title(sprintf(mytitle,t,E,P));
    drawnow;
end
end