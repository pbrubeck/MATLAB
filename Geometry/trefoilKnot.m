function S = trefoilKnot(n, m)
u=2*pi*(0:n-1)'/n;
C=[sin(u)+2*sin(2*u), cos(u)-2*cos(2*u), -sin(3*u)];
S=tubular(C, 0.5, m);
xx=S(:,:,1);
yy=S(:,:,2);
zz=S(:,:,3);
xx=[xx, xx(:,1)]; xx=[xx; xx(1,:)];
yy=[yy, yy(:,1)]; yy=[yy; yy(1,:)];
zz=[zz, zz(:,1)]; zz=[zz; zz(1,:)];
surf(xx,yy,zz);
camlight; shading interp; axis equal;
end