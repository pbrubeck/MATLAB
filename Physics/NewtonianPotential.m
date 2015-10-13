function G = NewtonianPotential(u, a, b)
% Solves Laplace equation for the scalar potential given the field's density 
% Assumes periodic boundary conditions
N=size(u);
P=b-a;
th=2i*pi/P;
gx=a+P/N(1)*(0:N(1)-1);
gy=a+P/N(2)*(0:N(2)-1);
gz=a+P/N(3)*(0:N(3)-1);
[xx,yy,zz]=meshgrid(gx, gy, gz);

i=th*[0:N(1)/2, -N(1)/2+1:-1];
j=th*[0:N(2)/2, -N(2)/2+1:-1];
k=th*[0:N(3)/2, -N(3)/2+1:-1];
[i2,j2,k2]=meshgrid(i.^2, j.^2, k.^2);
D2=i2+j2+k2;

px=exp(-a*i);
py=exp(-a*j);
pz=exp(-a*k);
phase=reshape(kron(kron(px,py),pz), N(1:3));
op=phase./D2; op(1,1,1)=0;


if(ndims(u)==4)
    u_hat(:,:,:,1)=fftn(u(:,:,:,1));
    u_hat(:,:,:,2)=fftn(u(:,:,:,2));
    u_hat(:,:,:,3)=fftn(u(:,:,:,3));
    G(:,:,:,1)=real(ifftn(op.*u_hat(:,:,:,1)));
    G(:,:,:,2)=real(ifftn(op.*u_hat(:,:,:,2)));
    G(:,:,:,3)=real(ifftn(op.*u_hat(:,:,:,3)));
    figure(1); axis equal;
    quiver3(xx,yy,zz,G(:,:,:,1),G(:,:,:,2),G(:,:,:,3));
else
    u_hat=fftn(u);
    G=real(ifftn(op.*u_hat));
    figure(1); axis equal
    M=a+(b-a)*(1/2-1./N);
    contourslice(xx,yy,zz,G,0,0,0,N(1));
    colormap(jet(128)); colorbar();
end
view(3);
end