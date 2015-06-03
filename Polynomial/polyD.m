function D = polyD(A, n)
% Returns the nth derivate of a set of polynomials.
if n==0
    D=A;
    return;
end
D=zeros(size(A,1), size(A,2)-n);
for i=1:size(A,2)-n
    k=1;
    for j=i:i+n-1
        k=k*j;
    end
    D(:,i)=k*A(:,i+n);
end
end