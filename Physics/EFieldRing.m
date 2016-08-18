function [] = EFieldRing(n, m)
R=1;
L=2*R;

x=linspace(-L, L, n);
t=2*pi/m*((0:m-1)+0.5);
[xx,zz]=meshgrid(x,x);
[xxx,zzz,phi]=meshgrid(x,x,t);

% Separation vector rr=r-r'
rr1=xxx-R*cos(phi);
rr2=-R*sin(phi);
rr3=zzz;
rr=rr1.^2+rr2.^2+rr3.^2;

% Integration
E1=sum(rr1./(rr).^(3/2), 3);
E2=sum(rr2./(rr).^(3/2), 3);
E3=sum(rr3./(rr).^(3/2), 3);

% Intensity image
EE=sqrt(E1.^2+0*E2.^2+E3.^2);
V=log2(1+mat2gray(EE));

id=(1+ceil(n/64):ceil(n/32):n-ceil(n/64));

figure(1); clf; hold on;
image([-L,L], [-L,L], 255*V);
set(gca, 'YDir', 'normal');
colormap(hot(256));
quiver(xx(id,id),zz(id,id),E1(id,id)./EE(id,id),E3(id,id)./EE(id,id),'w');
axis equal; axis off; hold off;


z=1;
[xx,yy]=meshgrid(x,x);
[xxx,yyy,phi]=meshgrid(x,x,t);

% Separation vector rr=r-r'
rr1=xxx-R*cos(phi);
rr2=yyy-R*sin(phi);
rr3=z;
rr=rr1.^2+rr2.^2+rr3.^2;

% Integration
E1=sum(rr1./(rr).^(3/2), 3);
E2=sum(rr2./(rr).^(3/2), 3);
E3=sum(rr3./(rr).^(3/2), 3);

% Intensity image
EE=sqrt(E1.^2+E2.^2+0*E3.^2);
V=log2(1+mat2gray(EE));

figure(2); clf; hold on;
image([-L,L], [-L,L], 255*V);
set(gca, 'YDir', 'normal');
colormap(hot(256));
quiver(xx(id,id),yy(id,id),E1(id,id)./EE(id,id),E2(id,id)./EE(id,id),'w');
axis equal; axis off; hold off;

end