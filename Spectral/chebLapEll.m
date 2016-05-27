function [F, G, Duu, Dvv, u, v] = chebLapEll(N,M)
% Matrices to calculate the laplacian on the unit disk
% usage: F*Lap(U)+Lap(U)*G=Duu*U+U*Dvv'
% domain r:(0,1], phi:(0,2*pi]
[Dvv,v]=periodicD2(M);
[Du,u]=chebD(2*N);
Duu=Du^2;
Duu=Duu(1:N,1:N)+Duu(1:N,end:-1:N+1);
u=u(1:N);
F=diag(cosh(2*u)/2);
G=diag(-cos(2*v)/2);
end
