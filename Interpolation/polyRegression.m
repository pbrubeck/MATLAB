function p = polyRegression(x, y, deg)
% Returns the polynomial that best fits the data.
A=ones(length(x), deg+1);
for j=1:deg
    A(:,j+1)=A(:,j).*x(:);
end
p=((A'*A)\A')*y;
end