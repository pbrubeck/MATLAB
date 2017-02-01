function S = torusKnot(p, q, n, m)
u=ainit(2*pi*(0:n-1)'/n, 3);
r=cos(q*u)+2;
C=[r*cos(p*u); r*sin(p*u); -sin(q*u)];
S=tubular(C,0.5,m);

xx=S([1:end,1],[1:end,1],1);
yy=S([1:end,1],[1:end,1],2);
zz=S([1:end,1],[1:end,1],3);
surf(xx,yy,zz);
camlight; shading interp; axis equal; view(2);
end