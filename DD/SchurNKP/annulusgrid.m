function [z0,quads,curv] = annulusgrid(R1,R2)
z0=kron([R1;R2],exp(1i*(0:3)'/4*2*pi));
quads=[1,5,2,6; 2,6,3,7; 3,7,4,8; 4,8,1,5];
curv=zeros(size(quads)); curv(:,1)=-1; curv(:,2)=-2;
end

