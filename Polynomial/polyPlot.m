function []=polyPlot(P, x, y)
if(nargin==2)
    y=Horner(P,x);
else
    y=Horner(P,y);
end
plot(x,y);
end