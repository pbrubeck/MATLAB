function [] = Young(m,n)
k=2*pi;
R=1;
b=4;
z=1000;

x=linspace(-z,z,n);
[xx,yy]=ndgrid(x);

E=zeros(n);
for i=0:m-1
    x0=b*cos(2*pi*i/m);
    y0=b*sin(2*pi*i/m);
    r=sqrt((xx-x0).^2+(yy-y0).^2+z^2);
    rho=k*R*hypot(xx-x0,yy-y0)./r;
    E=E+exp(1i*k*r).*jinc(rho);
end
I=real(conj(E).*E);

figure(1);
imagesc(x/z,x/z,I.');
colormap(hot(256));
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'YDir','normal','TickLabelInterpreter','latex','fontsize',14);
axis equal;
xlim([-1,1]); 
ylim([-1,1]);
xlabel('$x/z_R$');
ylabel('$y/z_R$');
%print('-depsc', sprintf('Young%02d',m));
end