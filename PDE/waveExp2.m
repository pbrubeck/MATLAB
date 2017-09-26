function [] = waveExp2(m,n)
if nargin<2
    n=m;
end

[Dx,x]=chebD(m);
[Dy,y]=chebD(n); y=y';
[xx,yy]=ndgrid(x,y);

E1=eye(m);
E2=eye(n);
Dxx=Dx*Dx;
Dyy=Dy*Dy;

opA=@(uu) Dxx*uu+uu*Dyy';

a=[1,1;1,1];
b=[0,0;0,0];
rd1=[1,m];
rd2=[1,n];
B1=diag(a(1,:))*E1(rd1,:)+diag(b(1,:))*Dx(rd1,:);
B2=diag(a(2,:))*E2(rd2,:)+diag(b(2,:))*Dy(rd2,:);

[gf,ps,kd,gb,dL,hyp] = elliptic(Dxx,E1,E2,Dyy,B1,B2,rd1,rd2);

f=@(x,y) real(sin((xx+1i*yy).^4))+exp(-100/2*(xx.^2+yy.^2));
u0=f(xx,yy);
w0=gf(kd(opA(u0)));
psi=u0-w0;
[U,Ut]=hyp(w0,-0.5*Dx*w0-0.5*w0*Dy');

xq=linspace(-1,1,m);
yq=linspace(-1,1,n);
[xxx,yyy]=ndgrid(xq, yq);
uuu=interpcheb(interpcheb(u0,xq,1),yq,2);

% Energy integral
[~,w1]=ClenshawCurtis(-1,1,n);
[~,w2]=ClenshawCurtis(-1,1,m);
E=@(uu,uut) w1*(uut.^2+(Dx*uu).^2+(uu*Dy').^2)*w2';

figure(1);
h1=surf(xxx,yyy,uuu);
colormap(jet(256)); colorbar;
shading interp; camlight;
axis manual;

T=30;
dt=0.05;
for t=0:dt:T
    uu=psi+U(t);
    
    disp(E(uu,Ut(t)));
    uuu=interpcheb(interpcheb(uu,xq,1),yq,2);
    set(h1,'ZData',uuu);
    drawnow;
end

end

