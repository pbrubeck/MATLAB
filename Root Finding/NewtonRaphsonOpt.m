function x = NewtonRaphsonOpt(f, x)
% Solves a system of non-linear equations
g=agrad(f,length(x));
H=ahess(f,length(x));
i=0;
y=1;
while(norm(y)>2*eps && i<60)
    vars=num2cell(x);
    y=g(vars{:});
    x=x-H(vars{:})\y;
    i=i+1;
end
end