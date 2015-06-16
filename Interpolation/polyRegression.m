function p = polyRegression(x, y, n)
% Returns the polynomial that best fits the data.
A=ones(length(x), n+1);
for j=1:n
    A(:,j+1)=A(:,j).*x(:);
end
p=(A'*A)\(A'*y(:));
end