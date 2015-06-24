function J = quadrature2d(f, x, y, t, w)
% Integrates f over the 2D region bounded by the functions x, y(x)
% using a standard quadrature rule for the [-1, 1] interval.
tic
n=length(t); n2=n*n;
xi=(x(2)+x(1))/2;
dx=(x(2)-x(1))/2;
xi=xi+dx*t;

lim=y(xi);
yj=(lim(2,:)+lim(1,:))/2;
dy=(lim(2,:)-lim(1,:))/2;

if(n>256)
    xi=gpuArray(xi);
    yj=gpuArray(yj);
    dy=gpuArray(dy);
end

xi=repeat(xi, n2);
yj=repeat(yj, n2);
dy=repeat(dy, n);
yj=yj+kron(dy,t);

dA=kron(dx*dy.*w, w);
J=gather(f(xi, yj)*dA(:));
toc
end

function X = repeat(x, N)
X=repmat(x, N/length(x), 1);
X=X(:).';
end