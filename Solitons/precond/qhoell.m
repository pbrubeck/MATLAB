function [] = qhoell(m,n)
a=20*sqrt(2);
b=15*sqrt(2);

omega=2/a;
lam=omega*0;

[xi,eta,jac,M,H,U,hshuff,J1,J2]=schrodell(m,n,a,b,lam);
c=sqrt(a^2-b^2);
xx=c*cosh(xi).*cos(eta);
yy=c*sinh(xi).*sin(eta);
ii=1:m;
jj=[1:n,1];

rr=hypot(yy,xx);
VH=(omega*rr).^2;
function uvu=expval(V,u)
    ju=J1*u*J2';
    jv=J1*V*J2';
    vu=jac.*jv.*ju;
    uvu=ju(:)'*vu(:);
end

nx=2;
ny=2;
u0 = HermitePsi([zeros(nx,1);1],xx*sqrt(omega)).*...
     HermitePsi([zeros(ny,1);1],yy*sqrt(omega));
u=u0/sqrt(M(u0,u0));

t=0;
setlatex();
mytitle='$z = %f$, $E/\\omega = %f$, $P = %f$';
figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-a,a]/sqrt(2));
ylim([-b,b]/sqrt(2));
pbaspect([a,b,sqrt((a^2+b^2)/2)]);
shading interp;
colormap(magma(256));
colorbar();
view(2);
E=real(H(u,u)+expval(VH,u))/(2*omega);
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
    E=real(H(u,u)+expval(VH,u))/(2*omega);
    P=real(M(u,u));
    title(sprintf(mytitle,t,E,P));
    drawnow;
end
end