function [y] = ShermanMorrison(Ainv, u, v, x)
z=Ainv(x);
w=Ainv(u);
y=z-(v'*z)/(1+v'*w)*w;
end

