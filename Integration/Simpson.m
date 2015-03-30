function j = Simpson(f, a, b, n)
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
j=3*h/8*(f(a)+sum+f(b));
end
