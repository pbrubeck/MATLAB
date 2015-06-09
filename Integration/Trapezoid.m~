function J = Trapezoid(f, a, b)
% Iterative trapezoidal method.
n=1;
h=b-a;
J=h*(f(a)+f(b))/2;
err=1;
while(err>1E-5)
    s=0;
    n=2*n;
    xi=a+h/2;
    for j=1:2:n-1
        s=s+f(xi);
        xi=xi+h;
    end
    h=h/2;
    aux=J/2+h*s;
    err=abs((J-aux)/aux);
    J=aux;
end
end