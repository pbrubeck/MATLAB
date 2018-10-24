function [] = incedemo(N,p,m1,m2,c,omega,L)

q=omega*c^2;
[x,w]=gaulob(-L,L,N);
[xx,yy]=ndgrid(x,x);
jac=w(:)*w(:)';

rr=hypot(yy,xx);
th=atan2(yy,xx);
zz=acosh((xx+1i*yy)/c);
xi=real(zz);
eta=imag(zz);

function P=mass(u,v)
    bv=jac.*v;
    P=u(:)'*bv(:);
end

uu=igbeam(xi,eta,rr,p,m1,m2,q,omega,@mass);

figure(1);
surf(xx,yy,real(uu));
shading interp;
colormap(magma(256));
colorbar();
view(2);

figure(2);
surf(xx,yy,imag(uu));
shading interp;
colormap(magma(256));
colorbar();
view(2);

figure(3);
surf(xx,yy,abs(uu).^2);
shading interp;
colormap(magma(256));
colorbar();
view(2);


figure(4);
surf(xx,yy,angle(uu));
shading interp;
colormap(hsv(256));
colorbar();
view(2);




end

