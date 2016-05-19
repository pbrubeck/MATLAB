function J = Trapezoid(f, a, b, n)
% Simple trapezoidal method.
h=(b-a)/n;
xi=a;
J=(f(a)+f(b))/2;
for j=1:1:n-1
   xi=xi+h;
   J=J+f(xi);
end
J=h*J;
end