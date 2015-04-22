function y = RungeKutta(f, yi, ti, h, n)
% Solves an initial value problem of the form dy/dt=f(t, y).
y=zeros(n);
for ii=1:n
    k1=f(ti, yi);
    k2=f(ti+h/2, yi+h/2*k1);
    k3=f(ti+h/2, yi+h/2*k2);
    k4=f(ti+h, yi+h*k3);
    yi=yi+h/6*(k1+2*k2+2*k3+k4);
    ti=ti+h;
    y(ii)=yi;
end
end