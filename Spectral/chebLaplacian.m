function [L, xx, yy] = chebLaplacian( N )
% Returns Laplacian matrix and tensor product grid (excludes the boundary)
[D, x]=chebD(N);
D2=D^2; D2=D2(2:N-1, 2:N-1);
I=eye(N-2);
L=kron(I,D2)+kron(D2,I);
[xx,yy]=meshgrid(x(2:N-1));
end