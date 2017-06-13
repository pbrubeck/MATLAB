function [f,df] = mapline(a,b)
% Map from [0,1] to a line segment connecting points a and b
f=@(t) (b-a)*t+a;
df=@(t) (b-a);
end