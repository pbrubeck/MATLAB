function [p, w] = NewtonPoly(x, y, t)
n=length(x);
w=y;
for i=1:n
    for j=1:i-1
        w(i)=(w(i)-w(j))/(x(i)-x(j));
    end
end

m=length(t);
p=zeros(m,1);
for i=1:m
    aux=0;
    for j=n:-1:2
        aux=(aux+w(j))*(t(i)-x(j-1));
    end
    p(i)=aux+w(1);
end

end