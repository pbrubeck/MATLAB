function [] = phasePlot(f, x0, x1, y0, y1, n)
x=linspace(x0, x1, n);
y=linspace(y0, y1, n);
[re, im]=meshgrid(x, y);
z=re+1i*im;

C=log(f(z(1:end)));
C=reshape(C,[n,n]);

map=hsv(256);
A=ind2rgb(uint8(128*(1+imag(C)/pi)), map);
B=mat2gray(real(C));
B=cat(3,B,B,B);

colormap(map);
image([x0,x1], [y0,y1], A.*B);

xlabel('Re(z)');
ylabel('Im(z)');
caxis([-pi,pi]);
colorbar('YTick', linspace(-pi, pi,5), ...
    'YTickLabel', {'-\pi','-\pi/2','0','\pi/2','\pi'});
end
