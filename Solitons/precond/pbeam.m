function [] = pbeam(m,n)
L=10;
u0=vnlsetp2(m,n,L);


beta=-1;
lam=0.5;
pot=@(r) -30*(besselj(1,5*r)).^2;
[U,H,P,Q,rr,th]=pnlse2(m,n,sqrt(2)*L,lam,pot);
xx=rr.*cos(th);
yy=rr.*sin(th);

u=u0;
t=0;
E=real(H(u,u)+Q(u,(beta/2)*abs(u).^2,u));
p=real(P(u,u));
display(E);

setlatex();
figure(1);
hs=surf(xx(:,[1:end,1]),yy(:,[1:end,1]),abs(u(:,[1:end,1])).^2);
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(sprintf('$z = %f$, $E = %f$, $P = %f$',t,E,p));
drawnow;


T=2*pi;
nframes=1000;
dt=T/nframes;
umax=zeros(nframes+1,1);
umax(1)=max(abs(u(:)).^2);
for i=1:nframes
    u=U(dt/2,u);
    u=u.*exp(-1i*beta*dt*(abs(u).^2));
    u=U(dt/2,u);
    t=t+dt;

    E=real(H(u,u)+Q(u,(beta/2)*abs(u).^2,u));
    p=real(P(u,u));

    set(hs,'ZData',abs(u(:,[1:end,1])).^2);
    title(sprintf('$z = %f$, $E = %f$, $P = %f$',t,E,p));
    drawnow;
    
    umax(i+1)=max(abs(u(:)).^2);
end

figure(3);
tt=linspace(0,T,nframes+1);
plot(tt,umax-mean(umax),'b');
title('$\max|\psi|-\langle\max|\psi|\rangle$');
xlabel('$z$');
end

