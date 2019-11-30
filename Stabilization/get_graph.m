function [iflux] = get_graph(itopo,vx,vy,wx,wy)
    tol=1e-6;
    nfaces=size(itopo,1);
    nel=size(itopo,2);    
    if(length(vx)==nel)
    vx=reshape(vx,1,1,nel);
    vy=reshape(vy,1,1,nel);
    end
    
    heat=zeros(nel,1);
    iflux=zeros(nfaces,nel);
    for e=1:nel
        ax=sum(wx(:,e)); 
        ay=sum(wy(:,e));
        iflux(1,e)=-wx(:,e)'*vx(:,1,e)/ax;
        iflux(2,e)=wx(:,e)'*vx(:,end,e)/ax;
        iflux(3,e)=-vy(1,:,e)*wy(:,e)/ay;
        iflux(4,e)=vy(end,:,e)*wy(:,e)/ay;
        
        vol=wx(:,e)'*ones(size(vx(:,:,e)))*wy(:,e);
        ux=wx(:,e)'*vx(:,:,e)*wy(:,e)/vol;
        uy=wx(:,e)'*vy(:,:,e)*wy(:,e)/vol;
        heat(e)=ux*ux+uy*uy;
    end

    %heat=-heat;
    %heat(:)=0;
    iflux=sign(iflux).*(abs(iflux)>tol);
    iflux(itopo==repmat(1:nel,nfaces,1))=0;
    
    igrad=zeros(nfaces,nel);
    for e=1:nel
    for f=1:nfaces
        ee=itopo(f,e);
        if(heat(e)-heat(ee)>tol)
            igrad(f,e)=1;
        elseif(heat(ee)-heat(e)>tol)
            igrad(f,e)=-1;
        end
    end
    end
    
    for e=1:nel
    for f=1:nfaces
        if(itopo(f,e)==e)
            iflux(f,e)=0;
        elseif(iflux(f,e)==0)
            iflux(f,e)=igrad(f,e);
        end
    end
    end
end

