function [] = DirichletBipolar(R1, R2, L, N)
N(1:2)=N;
a=sqrt((R1^2-R2^2)^2-8*L^2*(R1^2+R2^2)+16*L^4)/(4*L);
x1=-asinh(a/R1);
x2= asinh(a/R2);

[Dx,x]=chebD(N(1)); x=x1+(x2-x1)/2*(x+1); Dx=2/(x2-x1)*Dx; Dxx=Dx*Dx; 
[Dyy,y]=fourD2(N(2));
[xx,yy]=ndgrid(x,y);
zz=xx+1i*yy;
ww=a*coth(zz/2)-a*coth(x2)+L;
uu=real(ww);
vv=imag(ww);
J=abs(-a/2*csch(zz/2).^2).^2;
F=zeros(N);

bc=[1+0*y; -1+0*y];
RHS=J.*F-Dxx(:,[1,end])*bc;
psi=zeros(N);
psi([1,end],:)=bc;
psi(2:end-1,:)=sylvester(Dxx(2:end-1,2:end-1), Dyy', RHS(2:end-1,:));

h=mesh(uu(:,[1:end,1]),vv(:,[1:end,1]),psi(:,[1:end,1]));
colormap(jet(256)); shading interp;
set(h,'FaceColor','none'); view(2);
c=4*(R1+R2+L);% xlim([-c,c]); ylim([-c,c]);
xrange=2*c;
yrange=2*c;
zrange=max(psi(:))-min(psi(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
end