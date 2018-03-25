function [z,quads,curv] = quadmeshrefine(z,quads,curv,adj,bnd)
edges=[adj(:,1:2); bnd(:,1:2)];
nedges=size(edges,1);
ewns=[1,3; 2,4; 1,2; 3,4];
zold=z;
zcntr=mean(z(quads), 2);
zedge=zeros(size(edges,1),1);
for r=1:size(edges,1)
    zedge(r)=mean(z(quads(edges(r,2), ewns(edges(r,1),:))), 1);
end

z=[zedge; zcntr; zold];

nzedge=numel(zedge);
nzcntr=numel(zcntr);
nq=size(quads,1);

alledges=[adj(:,1:2); bnd(:,1:2); adj(:,3:4)];
mid=zeros(nq,4);
for i=1:4
    mid(alledges(alledges(:,1)==i,2), i) = mod(find(alledges(:,1)==i)-1, nedges)+1;
end

quads=quads+nzedge+nzcntr;
qnew=zeros(4*nq,4);
k=(1:nq)';
qnew(1:4:end,:)=[quads(:,1), mid(:,3), mid(:,1), k+nzedge];
qnew(2:4:end,:)=[mid(:,3), quads(:,2), k+nzedge, mid(:,2)];
qnew(3:4:end,:)=[mid(:,1), k+nzedge, quads(k,3), mid(:,4)];
qnew(4:4:end,:)=[k+nzedge, mid(:,2), mid(:,4), quads(:,4)];

quads=qnew;
curv=zeros(size(quads)); 
curv(:)=inf;
end

