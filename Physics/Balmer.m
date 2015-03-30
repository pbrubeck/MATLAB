function B=Balmer(n)
B=zeros(n,n);
for i=1:n
    for j=1:n
       B(i,j)=1/i^2-1/j^2; 
    end
end
end