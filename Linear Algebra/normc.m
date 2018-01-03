function [B] = normc(A, G)
% Normalizes columns of A with respect to metric G
if(nargin==2)
    if numel(G)==length(G)
        B=bsxfun(@rdivide, A, sqrt(G(:)'*(real(A).^2+imag(A).^2)));
    else
        B=bsxfun(@rdivide, A, sqrt(diag(A'*G*A)).');
    end
else
    B=bsxfun(@rdivide, A, sqrt(sum(real(A).^2+imag(A).^2)));
end
end