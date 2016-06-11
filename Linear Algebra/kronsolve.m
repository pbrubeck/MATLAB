function [X] = kronsolve(A, F)
% Solves multidimensional Sylvester equation 
% (kron(I,I,A1)+kron(I,A2,I)+kron(A3,I,I))*vec(X)=vec(F)
dim=ndims(F);
d=cell(dim,1);
V=cell(dim,1);
U=cell(dim,1);
if(~iscell(A))
    B=cell(dim,1);
    [B{1:dim}]=deal(A);
    A=B;
end
for i=1:dim
    [V{i}, d{i}]=eig(A{i},'vector');
end
if(dim==3)
    [di,dj,dk]=ndgrid(d{1},d{2},d{3});
    D=di+dj+dk;
else
    [di,dj]=ndgrid(d{1},d{2});
    D=di+dj;    
end
for i=1:dim
    U{i}=inv(V{i});
end
G=gekv(U, F, dim);
Z=G./D;
X=gekv(V, Z, dim);
end

function [x] = gekv(A, x, dim)
for i=1:dim
    x=permute(reshape(A{i}*reshape(x, size(A{i},1), []), size(x)), [2:dim 1]);
end
end