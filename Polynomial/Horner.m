function y = Horner(P, x)
% Efficient polynomial evaluation method.
y=zeros(size(P,1), length(x));
for j=size(P,2):-1:1
    y=bsxfun(@plus, P(:,j), bsxfun(@times, x, y));
end
end