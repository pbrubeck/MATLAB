function [icolor] = color_greedy(itopo,elems)
nel=size(itopo,2);
if(nargin==1)
    elems=1:nel;
end
elems=reshape(elems,1,[]);
icolor=zeros(nel,1);
for e=elems
    k=1;
    while any(k==icolor(itopo(:,e)))
        k=k+1;
    end
    icolor(e)=k;
end
end