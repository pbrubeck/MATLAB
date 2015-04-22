function [p] = Lagrange(x, y, t)
% Evaluates the Lagrange interpolation polynomial.
n=length(x);
w=ones(size(x));
for i=1:n
    for j=1:n
        if j~=i
            w(i)=w(i)*(x(i)-x(j)); 
        end  
    end
end
g=1;
p=0;
for k=1:n
    d=t-x(k);
    p=p+y(k)./(w(k)*d);
    g=g.*d;
end
p=g.*p;
end