function A = GaussJordan(A)
n=size(A,1);
for i=1:n
    r=i;
    while(A(r,i)==0 && r<n)
        r=r+1;
    end
    if(i~=r)
        A([i r],:)=A([r i],:);
    end
    if(A(i:i)~=1)
        A(i,i+1:end)=A(i,i+1:end)/A(i,i);
        A(i,i)=1;
    end
    for k=1:n
        if(i~=k)
            A(k,i+1:end)=A(k,i+1:end)-A(k,i)*A(i,i+1:end);
            A(k,i)=0;
        end
    end
end
end