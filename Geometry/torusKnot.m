function [] = torusKnot(p, q, n, m)
u=2*pi*(0:n-1)'/n;
r=cos(q*u)+2;
C=[r.*cos(p*u), r.*sin(p*u), -sin(q*u)];
tubular(C, 0.5, m);
end