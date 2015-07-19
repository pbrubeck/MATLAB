function x = fixedPoint(f, x0)
% Solves a system of non-linear equations of the kind x=f(x).
i=0;
err=1;
while(err>eps && i<100)
    x=f(x0);
    err=max(abs((x-x0)./x));
    x0=x;
    i=i+1;
end  
end