function [iflux] = get_graph(itopo,vx,vy,wx,wy)
    nfaces=size(itopo,1);
    nel=size(itopo,2);    
    if(length(vx)==nel)
    vx=reshape(vx,1,1,nel);
    vy=reshape(vy,1,1,nel);
    end
    
    iflux=zeros(nfaces,nel);
    for e=1:nel
        ax=sum(wx(:,e)); 
        ay=sum(wy(:,e));
        iflux(1,e)=-wx(:,e)'*vx(:,1,e)/ax;
        iflux(2,e)=wx(:,e)'*vx(:,end,e)/ax;
        iflux(3,e)=-vy(1,:,e)*wy(:,e)/ay;
        iflux(4,e)=vy(end,:,e)*wy(:,e)/ay;
    end
    iflux=sign(iflux);
    
    iheat=zeros(nel);
    iheat(floor(nel/2))=1;
    
    itemp=-2*ones(nfaces,nel);
    isweep=zeros(nel,1);
    src=find(iheat>0);
    q=length(src);
    isweep(1:q)=src;
    p1=1;
    p2=q;
    qlag=0;
    while (q>qlag && q<=nel)
        qlag=q;
        for p=p1:p2
            e=isweep(p);
            iheat(e)=2;
            for f=1:nfaces
                ee=itopo(f,e);
                if(iheat(ee)~=2)
                    itemp(f,e)=2;
                    if(iheat(ee)==0)
                        q=q+1;
                        isweep(q)=ee;
                        iheat(ee)=1;
                    end
                end
            end
        end
        p1=p2+1;
        p2=q;
    end
    for e=1:nel
        for f=1:nfaces
            if(itopo(f,e)==e)
                iflux(f,e)=0;
            elseif(iflux(f,e)==0)
                iflux(f,e)=itemp(f,e);
            end
        end
    end
end

