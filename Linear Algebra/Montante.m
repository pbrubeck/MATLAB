function A = Montante(A)
n=size(A,1);
p=1;
for i=1:n
    r=i;
    while(A(r,i)==0 && r<n)
        r=r+1;
    end
    if(r~=i)
        A([i r],:)=A([r i],:);
    end
    for k=1:n
        if(k~=i)
            A(k,:)=(A(i,i)*A(k,:)-A(k,i)*A(i,:))/p;
        end
    end
    p=A(i,i);
end  
A=A/p;
end