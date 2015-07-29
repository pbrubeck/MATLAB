function [] = NewtonianPotential(a, b, n)
%Solves Laplace equation for the scalar potential given the field's density 

P=b-a;
th=2i*pi/P;
gv=a+P/n*(0:n-1);
[xx,yy,zz]=meshgrid(gv);
i=th*[0:n/2, -n/2+1:-1];
j=th*[0:n/2, -n/2+1:-1];
k=th*[0:n/2, -n/2+1:-1];
[i2,j2,k2]=meshgrid(i.^2, j.^2, k.^2);
D2=(i2+j2+k2);

px=exp(-a*i);
py=exp(-a*j);
pz=exp(-a*k);
phase=reshape(kron(kron(px,py),pz),[n,n,n]);

u_hat=phase./D2;
u_hat(1,1)=0;
uu=real(ifftn(u_hat));
umax=max(uu(:));
umin=min(uu(:));

isovalues = umin+(umax-umin)*(0:0.01:1);
contourslice(xx,yy,zz,uu,0,0,0,isovalues);
view(3);
colormap(jet(128));
end

