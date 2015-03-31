function x = NewtonRaphson(f, J, x0)
% Solves a system of non-linear f(x)=0 equations given the Jacobian matrix.
ii=0;
err=1;
while(err>=1E-15 && ii<30)
    x=x0-J(x0)\f(x0);
    err=max(abs((x-x0)./x));
    x0=x;
    ii=ii+1;
    disp(transpose(x));
end
end