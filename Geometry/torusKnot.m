function S = torusKnot(p, q, n, m)
u=2*pi*(0:n-1)'/n;
r=cos(q*u)+2;
C=[r.*cos(p*u), r.*sin(p*u), -sin(q*u)];
S=tubular(C, 0.5, m);
xx=S(:,:,1);
yy=S(:,:,2);
zz=S(:,:,3);
xx=[xx, xx(:,1)]; xx=[xx; xx(1,:)];
yy=[yy, yy(:,1)]; yy=[yy; yy(1,:)];
zz=[zz, zz(:,1)]; zz=[zz; zz(1,:)];
surfl(xx,yy,zz,'light');
shading interp; axis equal;
end