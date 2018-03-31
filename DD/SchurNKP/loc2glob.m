function [GID,nid] = loc2glob(m, corners, adj)
% Local to global indexing
GID=zeros(m,m,size(corners,1));
rd=[1,m];
se=m-2;

% Corners
b=(corners==0);
corners=corners+se*size(adj,1);
nid=max(corners(:));
corners(b)=0;

for c=1:size(corners,1)
    GID([1,end],[1,end],c)=reshape(corners(c,:),[2,2]);
end

% Edges
for e=1:size(adj,1)
    eid=se*(e-1)+1:se*e;
    
    % Quad 1
    s=adj(e,1);
    q=adj(e,2);
    if(s>2)
        % Vertical edge
        ns=rd(s-2);
        dr=1-2*(GID(1,ns,q)>GID(m,ns,q));
        GID(2:end-1,ns,q)=eid(abs(min(dr,dr*se)):dr:abs(max(dr,dr*se)));
    else
        % Horizontal edge
        ew=rd(s);
        dr=1-2*(GID(ew,1,q)>GID(ew,m,q));
        GID(ew,2:end-1,q)=eid(abs(min(dr,dr*se)):dr:abs(max(dr,dr*se)));
    end
    
    % Quad 2
    s=adj(e,3);
    q=adj(e,4);
    if(s>2)
        % Vertical edge
        ns=rd(s-2);
        dr=1-2*(GID(1,ns,q)>GID(m,ns,q));
        GID(2:end-1,ns,q)=eid(abs(min(dr,dr*se)):dr:abs(max(dr,dr*se)));
    else
        % Horizontal edge
        ew=rd(s);
        dr=1-2*(GID(ew,1,q)>GID(ew,m,q));
        GID(ew,2:end-1,q)=eid(abs(min(dr,dr*se)):dr:abs(max(dr,dr*se)));
    end

end

end