function [] = isPerfect(n)
% Determines if the given number is a perfect number.

sum=1;
for i=2:floor(sqrt(n))
    if(mod(n,i)==0)
        sum=sum+i+n/i;
    end
end
if(sum==n && n~=1)
    disp([num2str(n), ' is a perfect number.']);
end

