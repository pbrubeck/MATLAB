function x = NewtonRaphson(f, x)
% Solves a system of non-linear equations
J=ajac(f,length(x));
i=0;
y=1;
while(norm(y)>2*eps && i<30)
    vars=num2cell(x);
    y=f(vars{:});
    x=x-J(vars{:})\y;
    i=i+1;
end
end