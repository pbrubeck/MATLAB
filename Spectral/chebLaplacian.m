function [L, xx, yy] = chebLaplacian( N )
% Returns Laplacian matrix and tensor product grid
[D, x]=chebD(N);
D2=D^2; D2=D2(2:N, 2:N);
I=eye(N-1);
L=kron(I,D2)+kron(D2,I);
[xx,yy]=meshgrid(x);
end