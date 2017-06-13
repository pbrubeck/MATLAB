function [f,df] = mapellipse(z,p,q)
% Map from [0,1] to an ellipse segment connecting points p and q
p=p-z;
q=q-z;

A=[real(p), -imag(p); imag(p), real(p)];
x=A\[real(q); imag(q)];
a=acos(x(1));
r=1i*p*x(2)/sqrt(1-x(1)^2);

f=@(t) z+p*cos(a*t)+r*sin(a*t);
df=@(t) a*(-p*sin(a*t)+r*cos(a*t));
end