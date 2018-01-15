n=256;
N=n+32;
xxx=gaulob(-1,1,N);

[D,x]=chebD(n);
[xx,yy]=ndgrid(x);

f=@(x,y) real(sin(10000*pi*(1./(4+x.*y)).^5));
uu=f(xx,yy);

E=interpchebfun2(n,xxx);

for i=1:10
    tic;
    uuu=E(E(uu,'notransp').','notransp').';
    toc
end

[xxx,yyy]=ndgrid(xxx);

figure(1);
surf(xx,yy,uu);
colormap(jet(256));
shading interp;
camlight;

figure(2);
surf(xxx,yyy,uuu);
colormap(jet(256));
shading interp;
camlight;

figure(3);
surf(xxx,yyy,uuu-f(xxx,yyy));
colormap(jet(256));
shading interp;
camlight;
