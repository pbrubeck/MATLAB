function [] = phasePlot(f, x0, x1, y0, y1, n)
x=linspace(x0, x1, n);
y=linspace(y1, y0, n);
[xx, yy]=meshgrid(x, y);
C=f(xx+1i*yy);

map=hsv(256);
hue=angle(C)/(2*pi);
hue=hue-floor(hue);
H=ind2rgb(uint8(255*hue), map);

G=abs(C);
gamma=1/log(max(G(:)));

p=0.1;
V=p+(1-p)*log2(1+mat2gray(G.^gamma));
V=cat(3, V, V, V);

image([x0,x1], [y0,y1], H.*V);
set(gca, 'YDir', 'normal');

xlabel('Re(z)');
ylabel('Im(z)');

colormap(map);
caxis([-pi,pi]);
colorbar('YTick', linspace(-pi, pi, 5), ...
    'YTickLabel', {'0','\pi/2','\pi','3/2\pi','2\pi'});
end
