function u = heatPDE(f, srcx, srcy, src, c, n)
% Solves the heat equation for an initial temperature distribution
t=linspace(-1, 1, n);
[x,y]=meshgrid(t, t);
x=gpuArray(x);
y=gpuArray(y);
u=reshape(f(x(1:end), y(1:end)), n, n);
u(srcx,srcy)=src;

filter=[0,1,0;1,-4,1;0,1,0];

map=jet(128);
h=imagesc(u);
set(gca, 'Units', 'pixels');
set(gca, 'YDir', 'normal');
set(gca, 'xlimmode','manual',...
         'ylimmode','manual',...
         'zlimmode','manual',...
         'climmode','manual',...
         'alimmode','manual');
set(gcf,'doublebuffer','off');

colormap(map);
nframes=10000;
for i=1:nframes
    A=conv2(u, filter);
    u=u+c*A(2:end-1, 2:end-1);
    u([1 end],:)=0;
    u(:,[1 end])=0;
    u(srcx,srcy)=src;
    if(mod(i,10)==1)
        set(h, 'CData', gather(u));
        drawnow;
    end
end
end