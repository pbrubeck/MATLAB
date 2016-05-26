function [R2, Drr, Dff, r, phi] = chebLapPol(N,M)
% Matrices to calculate the laplacian on the unit disk
% usage: R2*Lap(U)=Drr*U+U*Dff'
% domain r:(0,1], phi:(0,2*pi]
[Dff,phi]=periodicD2(M);
[D,r]=chebD(2*N);
Drr=(diag(r)*D)^2;
Drr=Drr(1:N,1:N)+Drr(1:N,end:-1:N+1);
r=r(1:N);
R2=diag(r.^2);
end

