function [z0,quads,curv] = Lgrid()
z0=[0; 1; 1+1i; 1i; -1+1i; -1; -1-1i; -1i];
quads=[3,4,2,1; 4,5,1,6; 1,6,8,7];
curv=zeros(size(quads)); curv(:,:)=inf;
end