function [] = phasePlot(f, x0, x1, y0, y1, n)
x=linspace(x0, x1, n);
y=linspace(y1, y0, n);
[xx, yy]=meshgrid(x, y);
C=f(xx+1i*yy);
map=hsv(256);
image([x0,x1], [y0,y1], complex2rgb(C,0.1));
set(gca, 'YDir', 'normal');

xlabel('Re(z)');
ylabel('Im(z)');

colormap(map);
caxis([-pi,pi]);
colorbar('YTick', linspace(-pi, pi, 5), ...
    'YTickLabel', {'0','\pi/2','\pi','3/2\pi','2\pi'});
end
