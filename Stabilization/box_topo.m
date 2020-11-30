function [itopo] = box_topo(nex,ney,nez)
ndim=nargin;
if(ndim==2)
    nez=1;
end
nfaces=2*ndim;
nel=nex*ney*nez;
itopo=zeros(nfaces,nel);
eid=reshape(1:nel,nex,ney,nez);
tt=eid; tt(2:end  ,:,:)=eid(1:end-1,:,:); itopo(1,:)=tt(:);
tt=eid; tt(1:end-1,:,:)=eid(2:end  ,:,:); itopo(2,:)=tt(:);
if(ndim>=2)
tt=eid; tt(:,2:end  ,:)=eid(:,1:end-1,:); itopo(3,:)=tt(:);
tt=eid; tt(:,1:end-1,:)=eid(:,2:end  ,:); itopo(4,:)=tt(:);
end
if(ndim>=3)
tt=eid; tt(:,:,2:end  )=eid(:,:,1:end-1); itopo(5,:)=tt(:);
tt=eid; tt(:,:,1:end-1)=eid(:,:,2:end  ); itopo(6,:)=tt(:);   
end
end