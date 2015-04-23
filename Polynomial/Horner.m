 function y = Horner(A, x)
% Efficient polynomial evaluation method.
y=zeros(size(A,1), length(x));
for j=size(A,2):-1:1
    y=bsxfun(@plus, A(:,j), bsxfun(@times, x, y));
end
end