function [] = phasePlot(f, x0, x1, y0, y1, n)
x=linspace(x0, x1, n);
y=linspace(y0, y1, n);
[xx, yy]=meshgrid(x, y);
C=f(xx+1i*yy);

gamma=1/log(max(abs(C(:))));

map=hsv(256);
A=ind2rgb(uint8(128+128*angle(C)/pi), map);
B=mat2gray(abs(C).^gamma);
B=cat(3, B, B, B);

image([x0,x1], [y0,y1], A.*B);
set(gca, 'YDir', 'normal');

xlabel('Re(z)');
ylabel('Im(z)');

colormap(map);
caxis([-pi,pi]);
colorbar('YTick', linspace(-pi, pi, 5), ...
    'YTickLabel', {'-\pi','-\pi/2','0','\pi/2','\pi'});
end
