function S = torusKnot(p, q, n, m)
u=2*pi*(0:n-1)'/n;
r=cos(q*u)+2;
C=[r.*cos(p*u), r.*sin(p*u), -sin(q*u)];
S=tubular(C, 0.5, m);
surf(S(:,:,1), S(:,:,2), S(:,:,3));
axis equal;
end