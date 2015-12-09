function C = Chebyshev(n)
% Returns the coeficient matrix of the first n Chebyshev polynomials.
C=zeros(n,n);
C(1,1)=1;
if(n<2) 
    return; 
end
C(2,2)=1;
for i=3:n
    C(i,1)=-C(i-2,1);
end
for i=3:n
    for j=3:n
        C(i,j)=2*C(i-1,j-1);  
    end
end
for j=2:n
    for i=j+1:n
        C(i,j)=2*C(i-1,j-1)-C(i-2,j);  
    end
end
end