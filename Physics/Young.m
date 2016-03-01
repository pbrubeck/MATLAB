function [] = Young(n)
k=40*pi;
R=1e-3;
d=0.3;
z=1;

[xx,yy]=meshgrid(linspace(-1,1,n));

I=zeros(n);
figure(1);
h=imagesc(I);
colormap(hot(256));
colorbar();

Nframes=20;
for m=1:Nframes
    E=zeros(n);
    for i=0:m-1
        x0=d*cos(2*pi*i/m);
        y0=d*sin(2*pi*i/m);
        rho=k*R*hypot(xx-x0,yy-y0)./sqrt(xx.^2+yy.^2+z^2);
        E=E+exp(1i*k*sqrt((xx-x0).^2+(yy-y0).^2+z^2)).*jinc(rho);
    end
    I=real(conj(E).*E);
    set(h,'CData',I);
    drawnow;
    pause(0.1);
end
end