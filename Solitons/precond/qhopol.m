function [] = qhopol(m,n,L,omega)
lam=0;
VL=@(r) (omega*r).^2;
[rr,th,jac,M,H,U]=schrodpol(m,n,L,lam,VL);
xx=rr.*cos(th);
yy=rr.*sin(th);
c=10;
zz=acosh((xx+1i*yy)/c);
xi=real(zz);
eta=imag(zz);
ii=1:m;
jj=[1:n,1];

coord=3;
if(coord==1)
nx=3; ny=3;
u0 = HermitePsi([zeros(nx,1);1],xx*sqrt(omega)).*...
     HermitePsi([zeros(ny,1);1],yy*sqrt(omega));
elseif(coord==2)
nr=7; 
l=5;
w=omega*rr.^2;
c=[zeros(1,(nr-abs(l))/2),1];
u0 = LaguerreL(c,abs(l),w).*exp(-w/2).*...
     (w.^(abs(l)/2)).*exp(1i*l*th);
elseif(coord==3)
p=3; m1=3; m2=3;
q=omega*c^2;
u0 = igbeam(xi,eta,rr,p,m1,m2,q,omega,M);
end
u=u0/sqrt(M(u0,u0));

t=0;
setlatex();
mytitle='$z = %f$, $E/\\omega = %f$, $P = %f$';
figure(1);
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L,L]/sqrt(2));
ylim([-L,L]/sqrt(2));
axis square;
shading interp;
colormap(magma(256));
colorbar();
view(2);
E=real(H(u,u))/(2*omega);
P=real(M(u,u));
title(sprintf(mytitle,t,E,P));
drawnow;

display(E);
display(omega);

nframes=100;
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