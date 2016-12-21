function [A1, A2, B1, r, t] = chebLapPol(N,M)
% Matrices to calculate the laplacian on the unit disk
% A1*U + U*A2' = B1*Lap(U)
% domain r:(0,1], t:(0,2*pi]
[D,r]=chebD(2*N);
A1=diag(r)*D+diag(r.^2)*D^2;
A1=A1(1:N,1:N)+A1(1:N,end:-1:N+1);
r=r(1:N);
B1=diag(r.^2);
[A2,t]=fourD2(M);
end