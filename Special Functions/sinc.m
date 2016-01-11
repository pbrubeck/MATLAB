function [y] = sinc(x)
t=(x==0);
y=(sin(x)+t)./(x+t);
end

