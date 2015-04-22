function J = Simpson(f, a, b, n)
% Returns the definite integral of the function f by Simpson's 3/8 rule
h=(b-a)/n;
xi=a;
sum=0;
for i=1:n-1
    xi=xi+h;
    if(mod(i,3)==0)
        sum=sum+2*f(xi);
    else
        sum=sum+3*f(xi);
    end
end
J=3*h*(f(a)+sum+f(b))/8;
end
