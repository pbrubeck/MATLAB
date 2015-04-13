function [p] = Lagrange(x, y, t)
% Uses barycentric interpolation formula
n=length(x);
w=ones(n,1);

% Compute barycentric reciprocal weights
for i=1:n
    for j=1:n
        if j~=i
            w(i)=w(i)*(x(j)-x(i)); 
        end  
    end
end
g=1;
p=0;
for k=1:n
    if t==x(k)
        p=y(k);
        g=1;
        break;
    else
        p=p+y(k)/(w*(t-x(k)));
        g=g*(t-x(k));
    end
end
p=g*p;

end

