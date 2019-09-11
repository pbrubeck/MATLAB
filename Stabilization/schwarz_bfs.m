function [isweep,icolor,bcs] = schwarz_bfs(itopo,iflux)
ndim=2;
nfaces=2*ndim;
nel=size(itopo,2);

isweep=zeros(nel,1);
icolor=zeros(nel,1);
bcs=zeros(nfaces,nel);
bcs(:)=2; % Overlap
bcs(iflux<0)=1; % Dirichlet, inflow


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
    e=1;
    nex=sqrt(nel); ney=nel/nex;
    e=floor((nex/2)+1)+(ney-1)*nex;
    q=q+1;
    isweep(q)=e;
    icolor(e)=k;
end

p1=1;
p2=q;
qlag=0;
while (q>qlag && q<nel)
    k=k+1;
    qlag=q;
    for p=p1:p2
    e=isweep(p);   
    for f=1:nfaces
        ee=itopo(f,e);
        if(iflux(f,e)>0 && ee~=e)
            if(icolor(ee)==0)
                q = q+1;
                isweep(q)=ee;
                icolor(ee)=k;
            end
        end
    end
    end
    p1=p2+1;
    p2=q;
end

end