function err = normError(f, g, x, w)
% Returns the norm of the error given a quadrature
y=f(x)-g(x);
err=sqrt((y.*y)*w(:));
end

