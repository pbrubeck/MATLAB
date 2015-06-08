function [] = modPlot(f, x0, x1, y0, y1, n)
%Generates a complex plot.
x=linspace(x0, x1, n);
y=linspace(y0, y1, n);

[re, im]=meshgrid(x, y);
s=(re+1i*im);
z=abs(f(s(1:end)));
z=reshape(z,[n,n]);
surf(x,y,z,'EdgeColor','none','LineStyle','none');
shading interp;
end

