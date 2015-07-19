function J = NewtonCotes(f, a, b, n)
% Returns the definite integral of the function f by Simpson's 1/3 rule
h=(b-a)/n;
xi=a;
J=0;
for i=1:n-1
    xi=xi+h;
    if(mod(i,2)==0)
        J=J+2*f(xi);
    else
        J=J+4*f(xi);
    end
end
J=h*(f(a)+f(b)+J)/3;
end