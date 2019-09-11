function [itopo] = box_topo(nex,ney)
ndim=2;
nfaces=2*ndim;
nel=nex*ney;
eid=reshape(1:nel,nex,ney);
itopo=zeros(nfaces,nel);
tt=eid; tt(2:end  ,:)=eid(1:end-1,:); itopo(1,:)=tt(:);
tt=eid; tt(1:end-1,:)=eid(2:end  ,:); itopo(2,:)=tt(:);
tt=eid; tt(:,2:end  )=eid(:,1:end-1); itopo(3,:)=tt(:);
tt=eid; tt(:,1:end-1)=eid(:,2:end  ); itopo(4,:)=tt(:);
end