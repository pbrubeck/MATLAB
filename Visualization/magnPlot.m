function [] = magnPlot(f, x0, x1, y0, y1, n)
%Generates a complex plot.
x=linspace(x0, x1, n);
y=linspace(y0, y1, n);

[re, im]=meshgrid(x, y);

s=re+1i*im;
z=arrayfun(@(u) abs(f(u)), s);
surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
shading interp
end

