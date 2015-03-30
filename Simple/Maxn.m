function max=Maxn(arr,n)
% Maxn Returns the n largest elements from the given array.

max(1:n)=-realmax;
for i=1:length(arr)
    j=1;
    while j<=n
        if arr(i)>=max(j)
            for k=n:-1:j+1
                max(k)=max(k-1);
            end
            max(j)=arr(i);
            j=j+n;
        end
        j=j+1;
    end
end