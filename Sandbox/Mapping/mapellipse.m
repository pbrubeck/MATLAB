function [f,df] = mapellipse(z,p,q)
% Map from [0,1] to an ellipse segment connecting points p and q
p=p-z;
q=q-z;

a=pi/2;
f=@(t) z+p*cos(a*t)+q*sin(a*t);
df=@(t) a*(-p*sin(a*t)+q*cos(a*t));
end