function [pos] = ddpatches( topo )
% Computes the center of each patch of the domain decomposition
ndoms=size(topo,1);
EWNS=[1,-1,1i,-1i];
pos=zeros(ndoms,1);
dep=false(ndoms+1,1);
dep([1,end])=1;
k=1;
while(prod(dep)==0)
    m=min(topo(k,dep(topo(k,:))==0));
    pos(m)=pos(k)+2*EWNS(topo(k,:)==m);
    dep(m)=1;
    pot=find(dep(1:ndoms));
    T=topo(dep(1:ndoms),:);
    k=min(pot(prod(reshape(dep(T(:)),size(T)),2)==0));
end
end