function [B] = normc(A, w)
% Normalizes columns of A with respect to weight w
if(nargin==2)
    B=bsxfun(@rdivide, A, sqrt(w(:)'*(real(A).^2+imag(A).^2)));
else
    B=bsxfun(@rdivide, A, sqrt(sum(real(A).^2+imag(A).^2)));
end
end