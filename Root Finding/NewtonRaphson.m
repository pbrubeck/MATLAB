function x = NewtonRaphson(f, J, x)
% Solves a system of non-linear equations given the Jacobian matrix.
i=0;
y=1;
while(norm(y)>2*eps && i<30)
    y=f(x);
    x=x-J(x)\y;
    i=i+1;
end
end