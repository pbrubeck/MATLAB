function [B,a] = normc(A, G)
% Normalizes columns of A with respect to metric G
if(nargin==2)
    if numel(G)==length(G)
        a = sqrt(G(:)'*(real(A).^2+imag(A).^2));
        B=bsxfun(@rdivide, A, a);
    else
        a = sqrt(sum(A.*(G*A)),2).';
        B=bsxfun(@rdivide, A, a);
    end
else
    B=reshape(A,[],size(A,ndims(A)));
    a=sqrt(sum(real(B).^2+imag(B).^2));
    a(a==0)=1;
    B=bsxfun(@rdivide, B, a);
    B=reshape(B,size(A));
end
end