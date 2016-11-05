function psi=kmbreather(b,x,t)
% Kuznetsov-Ma breather, period 2*pi/b
a=(1+sqrt(1+b^2))/4; w=2*sqrt(2*a-1);
psi=(1+(2*(2*a-1)*cos(b*t)+1i*b*sin(b*t))./(sqrt(2*a)*cosh(w*x)+cos(b*t)));
end