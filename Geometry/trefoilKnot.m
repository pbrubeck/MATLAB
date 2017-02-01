function S = trefoilKnot(n, m)
u=ainit(2*pi*(0:n-1)'/n, 3);
C=[sin(u)+2*sin(2*u); cos(u)-2*cos(2*u); -sin(3*u)];
S=tubular(C,0.5,m);

xx=S([1:end,1],[1:end,1],1);
yy=S([1:end,1],[1:end,1],2);
zz=S([1:end,1],[1:end,1],3);
surf(xx,yy,zz);
camlight; shading interp; axis equal; view(2);
end