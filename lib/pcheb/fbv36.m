function f = fbv36(x,y)
%FBV36  Boundary-value test function (Trefethen, Problem 36).

    f = (x==1)       .* sin(3*pi*y)/5 ...
      + (y==1 & x<0) .* sin(pi*x).^4;

end

