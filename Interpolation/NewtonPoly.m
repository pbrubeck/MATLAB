function [p, a] = NewtonPoly(x, y, t)
% Evaluates the Newton interpolation polynomial.
n=length(x);
a=y;
for i=1:n
    for j=1:i-1
        a(i)=(a(i)-a(j))/(x(i)-x(j));
    end
end
p=zeros(size(t));
for j=n:-1:2
    p=(p+a(j)).*(t-x(j-1));
end
p=p+a(1);
end