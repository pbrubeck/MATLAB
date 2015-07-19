function p = leastSquares(x, y, deg)
% Returns the polynomial that best fits the data.
n=deg+1;
A=zeros(n,n+1);
xn=ones(size(x));
for i=1:2*n
    sx=sum(xn);
    if i<=n
        A(i, n+1)=xn*y(:);
    end
    row=min(i,n);
    col=1+i-row;
    for j=0:n-abs(i-n)-1
        A(row-j,col+j)=sx;
    end
    xn=xn.*x;
end
A=rref(A);
p=transpose(A(:,n+1));
end