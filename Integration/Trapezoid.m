function J = Trapezoid(f, a, b, n)
h=(b-a)/n;
xi=a;
J=(f(a)+f(b))/2;
for i=2:n
    xi=xi+h;
    J=J+f(xi);
end
J=J*h;
end

