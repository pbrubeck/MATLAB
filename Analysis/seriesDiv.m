function c=seriesDiv(a, b)
m=length(a);
n=length(b);
c=zeros(size(a));
for i=1:n
    t=a(i);
    for j=max(1,i-n+1):i
        t=t-b(i-j+1)*c(j);
    end
    c(i)=t/b(1);
end
end
