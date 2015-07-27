function [V, A] = surfVolume(S)
% Calculates the volume and surface area of a closed, smooth 2-manifold.
du=2*pi/size(S,1);
dv=2*pi/size(S,2);

Su=spectralD(S,1,1);
Sv=spectralD(S,1,2);

N=cross(Su, Sv);
dA=sqrt(sum(N.*conj(N), 3));
A=abs(sum(dA(:)))*du*dv;

uN=bsxfun(@ldivide, dA, N);
dV=dot(uN, S, 3);
V=abs(sum(dV(:)))*du*dv/3;
end