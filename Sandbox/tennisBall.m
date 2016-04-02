function [] = tennisBall(N)
figure(1);
[h,ph,th]=sphericalPlot(2*N,N);

A=0.5;
gamma=1/4;
t=cos(2*ph);
u=sign(t).*(A*abs(t).^gamma);
C=sign(2/pi*th-1-u);

alpha(0.6);
set(h,'CData',C);
shading interp;
colormap([1 0 0; 0 0 1])
end