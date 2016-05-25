function [V, A] = surfVolume(S)
% Calculates the volume and surface area of a closed, smooth 2-manifold.
du=2*pi/size(S,1);
dv=2*pi/size(S,2);

N=cross(fftD(S,1,1), fftD(S,2,1));
dA=sqrt(dot(N, N, 3));
A=sum(dA(:))*du*dv;

uN=bsxfun(@rdivide, N, dA);
dV=dot(uN, S, 3);
V=abs(sum(dV(:)))*du*dv/3;
end