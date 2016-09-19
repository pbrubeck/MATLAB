function [D] = complexD(z)
% Differentiation matrix on given nodes
N=length(z);
Z=repmat(z(:), [1, N]);
dZ=Z-Z.'+eye(N);
p=prod(dZ,2);
D=(p*(1./p).')./dZ-diag(sum(1./dZ));
end