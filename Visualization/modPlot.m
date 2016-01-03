function [] = modPlot(f, x0, x1, y0, y1, n)
%Generates a complex plot.
x=linspace(x0, x1, n);
y=linspace(y0, y1, n);

[re, im]=meshgrid(x, y);
s=re+1i*im;
z=real(log(f(s)));
surf(x,y,z,'EdgeColor','none','LineStyle','none');
colormap(jet(256));
xlabel('Re(z)');
ylabel('Im(z)');
shading interp;
end

