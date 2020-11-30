function x = NewtonRaphsonOpt(f, x)
% Solves a system of non-linear equations
g=agrad(f,length(x));
H=ahess(f,length(x));
i=0;
err=1;
tol=10*eps;
while(err>tol && i<60)
    vars=num2cell(x);
    y=g(vars{:});
    dx=H(vars{:})\y;
    err=norm(y);
    x=x-dx;
    i=i+1;
end
end