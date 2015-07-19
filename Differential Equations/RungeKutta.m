function y = RungeKutta(f, yi, ti, h, n)
% Solves an initial value problem of the form dy/dt=f(t, y).
y=zeros([size(yi) n]);
for i=1:n
    k1=h*f(ti, yi);
    k2=h*f(ti+h/2, yi+k1/2);
    k3=h*f(ti+h/2, yi+k2/2);
    k4=h*f(ti+h, yi+k3);
    yi=yi+(k1+2*k2+2*k3+k4)/6;
    ti=ti+h;
    y(i)=yi;
end
end