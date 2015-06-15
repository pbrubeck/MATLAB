function x = NewtonRaphson(f, J, x0)
% Solves a system of non-linear equations given the Jacobian matrix.
i=0;
err=1;
while(err>eps && i<30)
    x=x0-J(x0)\f(x0);
    err=max(abs((x-x0)./x));
    x0=x;
    i=i+1;
end
end
