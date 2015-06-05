function D = polyD(P, n)
% Returns the nth derivate of a set of polynomials.
if(nargin==1)
    n=1;
end
if n==0
    D=P;
    return;
end
D=zeros(size(P,1), size(P,2)-n);
for i=1:size(P,2)-n
    k=1;
    for j=i:i+n-1
        k=k*j;
    end
    D(:,i)=k*P(:,i+n);
end
end
