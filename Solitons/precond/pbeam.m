function [] = pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f)

L=max(sqrt((xx(:).^2+yy(:).^2)/2));
ii=1:size(u,1);
jj=[1:size(u,2),1];

function [ufu]=expval(f,u)
    ju=J1*u*J2';
    u2=abs(ju).^2;
    fu=jac.*f(u2).*ju;
    ufu=ju(:)'*fu(:);
end

t=0;
E=real(H(u,u)+expval(f,u))/2;
P=real(M(u,u));
display(E);

figure(1);
setlatex();
mytitle='$z = %f$, $E = %f$, $P = %f$';
h1=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
xlim([-L,L]);
ylim([-L,L]);
colormap(magma(256));
colorbar();
shading interp;
axis square;
view(2);
title(sprintf(mytitle,t,E,P));


figure(2);
h2=surf(xx(ii,jj),yy(ii,jj),angle(u(ii,jj)));
xlim([-L,L]);
ylim([-L,L]);
colormap(hsv(256));
colorbar();
shading interp;
axis square;
view(2);
title(sprintf(mytitle,t,E,P));

drawnow;

dt=T/nframes;
umax=zeros(nframes+1,1);
umax(1)=max(abs(u(:)).^2);
for i=1:nframes
    u=U(dt/2,u);
    u=u.*exp(-2i*dt*f(abs(u).^2));
    u=U(dt/2,u);
    t=t+dt;

    E=real(H(u,u)+expval(f,u))/2;
    P=real(M(u,u));
    set(h1,'ZData',abs(u(ii,jj)).^2);
    title(get(1,'CurrentAxes'),sprintf(mytitle,t,E,P));
    set(h2,'ZData',angle(u(ii,jj)));
    title(get(2,'CurrentAxes'),sprintf(mytitle,t,E,P));
    drawnow;
    
    umax(i+1)=max(abs(u(:)).^2);
end

figure(3);
tt=linspace(0,T,nframes+1);
plot(tt,umax-mean(umax),'b');
title('$\max|\psi|-\langle\max|\psi|\rangle$');
xlabel('$z$');
end