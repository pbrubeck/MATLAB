function y = EulerMethod(f, yi, ti, h, n)
% Solves an initial value problem of the form dy/dt=f(t, y).
y=zeros(n);
for i=1:n
    yi=yi+h*f(ti, yi);
    ti=ti+h;
    y(i)=yi;
end
end