function [] = phasePlot(f, x0, x1, y0, y1, n)
hmap(1:256,1) = linspace(0,1,256); 
hmap(:,[2 3]) = 1; %brightness 
huemap = hsv2rgb(hmap); 
colormap(huemap)

x=linspace(x0, x1, n);
y=linspace(y0, y1, n);
[re, im]=meshgrid(x, y);
C=re+1i*im;
image(length(huemap)/2*(angle(f(C))/pi+1));
colorbar;
end

