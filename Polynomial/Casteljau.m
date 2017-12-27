function [ y ] = Casteljau(a, x)
% Evaluation of Bernstein polynomial
n=length(a);
A=repmat(a(:)', numel(x), 1);
for j=2:n
    A(:,1:n-j+1)=A(:,1:n-j+1)+bsxfun(@times, x(:), diff(A(:,1:n-j+2),1,2));
end
y=reshape(A(:,1),size(x));
end

