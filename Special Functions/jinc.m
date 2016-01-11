function [j] = jinc(a,x)
t=(x==0);
j=(besselj(a,x)+t)./(x+t);
end

