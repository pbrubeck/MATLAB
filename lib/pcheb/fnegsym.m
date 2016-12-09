function f = fnegsym(x,y)
%FNEGSYM  Antisymmetric test function satisfying f(y,x) = f(-y,-x) = -f(x,y).

    f = x.^2 - y.^2;

end
