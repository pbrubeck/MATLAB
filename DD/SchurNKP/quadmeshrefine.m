function [z,quads,curv] = quadmeshrefine(z,quads,curv,adj,bnd)
edges=[adj(:,1:2); bnd(:,1:2)];
nedges=size(edges,1);
nquads=size(quads,1);

ewns=[1,3; 2,4; 1,2; 3,4];
zedge=zeros(size(edges,1),1);
for r=1:size(edges,1)
    pts=z(quads(edges(r,2), ewns(edges(r,1),:)));
    rad=curv(edges(r,2), edges(r,1));
    zedge(r)=mean(pts);
    if ~isinf(rad)
        d=diff(pts)/2;
        s=1i*sign(rad*d);
        h=sqrt(rad^2-abs(d)^2);
        zedge(r)=zedge(r)+s*(abs(rad)-h);
    end
end
zcntr=mean(z(quads), 2);
for q=1:nquads
    for j=1:4
        rad=curv(q, j);
        if ~isinf(rad)
            pts=z(quads(q, ewns(j,:)));
            d=diff(pts)/2;
            s=1i*sign(rad.*d);
            h=sqrt(rad.^2-abs(d).^2);
            zcntr(q)=zcntr(q)+s.*(abs(rad)-h)/2;
        end
    end
end

z=[zedge; zcntr; z];

nzedge=numel(zedge);
nzcntr=numel(zcntr);


alledges=[adj(:,1:2); bnd(:,1:2); adj(:,3:4)];
mid=zeros(nquads,4);
for i=1:4
    mid(alledges(alledges(:,1)==i,2),i) = 1+mod(find(alledges(:,1)==i)-1, nedges);
end

quads=quads+nzedge+nzcntr;
qnew=zeros(4*nquads,4);
k=(1:nquads)';
qnew(1:4:end,:)=[quads(:,1), mid(:,3), mid(:,1), k+nzedge];
qnew(2:4:end,:)=[mid(:,3), quads(:,2), k+nzedge, mid(:,2)];
qnew(3:4:end,:)=[mid(:,1), k+nzedge, quads(:,3), mid(:,4)];
qnew(4:4:end,:)=[k+nzedge, mid(:,2), mid(:,4), quads(:,4)];
quads=qnew;

curvnew=zeros(size(quads));
curvnew(1:4:end,[1,3])=curv(:,[1,3]);
curvnew(2:4:end,[2,3])=curv(:,[2,3]);
curvnew(3:4:end,[1,4])=curv(:,[1,4]);
curvnew(4:4:end,[2,4])=curv(:,[2,4]);

curvnew(1:4:end,2)=mean(curv(:,[1,2]),2);
curvnew(1:4:end,4)=mean(curv(:,[3,4]),2);
curvnew(2:4:end,1)=mean(curv(:,[1,2]),2);
curvnew(2:4:end,4)=mean(curv(:,[3,4]),2);
curvnew(3:4:end,2)=mean(curv(:,[1,2]),2);
curvnew(3:4:end,3)=mean(curv(:,[3,4]),2);
curvnew(4:4:end,1)=mean(curv(:,[1,2]),2);
curvnew(4:4:end,3)=mean(curv(:,[3,4]),2);

curvnew(curvnew==0)=inf;
curv=curvnew;
end