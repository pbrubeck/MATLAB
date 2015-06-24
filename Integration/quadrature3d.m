function J = quadrature3d(f, x, y, z, t, w)
% Integrates f over the 3D region bounded by the functions x, y(x), z(x,y)
% using a standard quadrature rule for the [-1, 1] interval.
n=length(t); n2=n*n; n3=n2*n;
xi=(x(2)+x(1))/2;
dx=(x(2)-x(1))/2;
xi=xi+dx*t;

lim=y(xi);
yj=(lim(2,:)+lim(1,:))/2;
dy=(lim(2,:)-lim(1,:))/2;

xi=repeat(xi, n2);
yj=repeat(yj, n2);
dy=repeat(dy, n);
yj=yj+kron(dy,t);

lim=z(xi, yj);
zk=(lim(2,:)+lim(1,:))/2;
dz=(lim(2,:)-lim(1,:))/2;

if(n>32)
    xi=gpuArray(xi);
    yj=gpuArray(yj);
    zk=gpuArray(zk);
    dz=gpuArray(dz);
end

xi=repeat(xi, n3);
yj=repeat(yj, n3);
zk=repeat(zk, n3);
dz=repeat(dz, n2);
zk=zk+kron(dz,t);

dA=kron(dx*dy.*w, w);
dV=kron(dA.*dz, w);

J=gather(f(xi, yj, zk)*dV(:));
end

function X = repeat(x, N)
X=repmat(x, N/length(x), 1);
X=X(:).';
end