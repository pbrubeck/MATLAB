function C=Chebyshev(n)
%Computes a 2n+1 by 2n+1 Chebyshev matrix

k=2*n+1;

C=zeros(k,k);
C(1,1)=1;
C(2,2)=1;
for i=3:k
    C(i,1)=-C(i-2,1);
end

for i=3:k
    for j=3:k
      C(i,j)=2*C(i-1,j-1);  
    end
end

for j=2:k
    for i=j+1:k
      C(i,j)=2*C(i-1,j-1)-C(i-2,j);  
    end
end

end

