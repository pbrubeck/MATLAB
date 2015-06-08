function [] = magnPlot(f, x0, x1, y0, y1, n)
%Generates a complex plot.
x=linspace(x0, x1, n);
y=linspace(y0, y1, n);

[re, im]=meshgrid(x, y);
s=re+1i*im;
C=(re+1i*im);
z=abs(f(reshape(C,[1,n*n])));
z=reshape(z,[n,n]);
surf(x,y,z,'EdgeColor','none','LineStyle','none');
shading interp;
end

