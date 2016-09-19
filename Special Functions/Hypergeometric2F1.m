function [ F ] = Hypergeometric2F1(a,b,c,z)
N=128;
if(isinteger(a) && a<0)
    N=-a+1;
end
if(isinteger(b) && b<0)
    N=-b+1;
end

F=zeros(size(z));
for k=N:-1:1
    F=1+(a+k-1)*(b+k-1)/((c+k-1)*k)*(z.*F);
end
end