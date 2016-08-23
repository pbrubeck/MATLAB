function [B] = normc(A, w)
% Normalizes columns of A with respect to weight w
B=bsxfun(@rdivide, A, sqrt(w(:)'*(real(A).^2+imag(A).^2)));
end