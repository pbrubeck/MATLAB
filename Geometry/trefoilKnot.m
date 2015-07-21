function S = trefoilKnot(n, m)
u=2*pi*(0:n-1)'/n;
C=[sin(u)+2*sin(2*u), cos(u)-2*cos(2*u), -sin(3*u)];
S=tubular(C, 0.5, m);
mesh(S(:,:,1), S(:,:,2), S(:,:,3));
axis equal;
end