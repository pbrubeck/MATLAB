function J = Simpson(f, a, b, n)
% Returns the definite integral of the function f by Simpson's 3/8 rule
h=(b-a)/n;
xi=a;
J=0;
for i=1:n-1
    xi=xi+h;
    if(mod(i,3)==0)
        J=J+2*f(xi);
    else
        J=J+3*f(xi);
    end
end
J=3*h*(f(a)+f(b)+J)/8;
end
