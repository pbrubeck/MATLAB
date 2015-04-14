function [p, a] = NewtonPoly(x, y, t)
% Evaluates the Newton interpolation polynomial.
n=length(x);
a=y;
for i=1:n
    for j=1:i-1
        a(i)=(a(i)-a(j))/(x(i)-x(j));
    end
end
m=length(t);
p=zeros(m,1);
for i=1:m
    aux=0;
    for j=n:-1:2
        aux=(aux+a(j))*(t(i)-x(j-1));
    end
    p(i)=aux+a(1);
end
end