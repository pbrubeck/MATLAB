function [isweep,icolor] = toposort_loops(itopo,iflux)
nfaces=size(itopo,1);
nel=size(itopo,2);

isweep=zeros(nel,1);
icolor=zeros(nel,1);
ifail=zeros(nel,1);

k=1;
q=0;
% Detect root
for e=1:nel
    if(all(iflux(:,e)>=0))
        q=q+1;
        isweep(q)=e;
        icolor(e)=k;
    end
end

if(q==0)
    % could not find a root node
    % define arbitrary root
    [fmax,e]=max(sum(min(0,iflux),1));
    q=q+1;
    isweep(q)=e;
    icolor(e)=k;
else
    % apply greedy coloring of root (only needed for DG)
    % FIXME greedy coloring of all neighbors
    iroot=find(icolor==1); 
    icolor=color_greedy(itopo,iroot);
    iroot=find(icolor==1);
    q=length(iroot);
    isweep(1:q)=iroot;
    icolor(icolor>1)=0;
end

t=0;
p1=1;
p2=q;
while(q<nel) 
    if(p2<p1) % frontier is empty
        % break loop closest to the root
        j=1;
        e=ifail(j);
        while(j<=t && icolor(e)>0)
            j=j+1;
            e=ifail(j);
        end
        k=-icolor(e);
        q=q+1;
        isweep(q)=e;
        icolor(e)=k;
        p1=q;
        p2=q;
        t=0;
    end
    k=k+1;
    for p=p1:p2
    e=isweep(p);   
    for f=1:nfaces
        ee=itopo(f,e);
        if(iflux(f,e)>=0 && ee~=e && ...
            all(icolor(itopo(:,ee))<k))
            ff=find(itopo(:,ee)==e);
            itopo(ff,ee)=ee;
            iflux(ff,ee)=0;
            itopo(f,e)=e;
            iflux(f,e)=0;
            % BFS with toposort priority
            if(icolor(ee)<=0 && all(iflux(:,ee)>=0))
                q = q+1;
                isweep(q)=ee;
                icolor(ee)=k;
            elseif(icolor(ee)==0)
                t = t+1;
                ifail(t)=ee;
                icolor(ee)=-k;
            end
        end
    end
    end
    p1=p2+1;
    p2=q;
end
[~,isweep]=sort(icolor);
end