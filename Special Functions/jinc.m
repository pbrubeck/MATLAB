function [j] = jinc(x)
t=(x==0);
j=(besselj(1,x)+t)./(x+t);
end

