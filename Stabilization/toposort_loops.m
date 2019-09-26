function [isweep,icolor] = toposort_loops(itopo,iflux)
ndim=2;
nfaces=2*ndim;
nel=size(itopo,2);

isweep=zeros(nel,1);
icolor=zeros(nel,1);
ifail=zeros(nel,1);

k=1;
q=0;
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
        if(iflux(f,e)>0 && ee~=e)
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



function [itopo,iflux]=cut_loops(itopo,iflux,i1,i2)
    if(i2>i1)
        % TODO check if DAG
        j1=floor((i1+i2-1)/2);
        j2=j1+1;
        nfaces=size(itopo,1);
        nab=0;
        for e=i1:j1
            for f=1:nfaces
                ee=itopo(f,e);
                if(j2<=ee && ee<=i2 && iflux(f,e)>0)
                    nab=nab+1;
                end
            end
        end
        nba=0;
        for e=j2:i2
            for f=1:nfaces
                ee=itopo(f,e);
                if(i1<=ee && ee<=j1 && iflux(f,e)>0)
                    nba=nba+1;
                end
            end
        end
        if(nba<nab)
            n1=j2; n2=i2;
            m1=i1; m2=j1;
        else    
            m1=j2; m2=i2;
            n1=i1; n2=j1;
        end
        for e=m1:m2
            for f=1:nfaces
                ee=itopo(f,e);
                if(n1<=ee && ee<=n2 && iflux(f,e)>0)
                    ff=(itopo(:,ee)==e);
                    iflux(ff,ee)=0;
                    itopo(ff,ee)=ee;
                    iflux(f,e)=0;
                    itopo(f,e)=e;
                end
            end
        end
        [itopo,iflux]=cut_loops(itopo,iflux,m1,m2);
        [itopo,iflux]=cut_loops(itopo,iflux,n1,n2);
    end
end

