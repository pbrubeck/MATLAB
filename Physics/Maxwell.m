function [] = Maxwell(rho, J)
% Solves Maxwell's (static) equations to find E and B using Fourier
% Spectral Methods
a=-pi; b=pi;
N=size(rho);
P=b-a;
th=2i*pi/P;
x=a+P/N(1)*(0:N(1)-1);
y=a+P/N(2)*(0:N(2)-1);
z=a+P/N(3)*(0:N(3)-1);
[xx,yy,zz]=ndgrid(x,y,z);

i=th*[0:N(1)/2, -N(1)/2+1:-1];
j=th*[0:N(2)/2, -N(2)/2+1:-1];
k=th*[0:N(3)/2, -N(3)/2+1:-1];
[ii,jj,kk]=ndgrid(i, j, k);
omega=cat(4, ii, jj, kk); % Frequency space
[i2,j2,k2]=ndgrid(i.^2, j.^2, k.^2);
D2=-(i2+j2+k2); % Laplacian operator

op=1./D2; op(1,1,1)=0; % Newtonian potential operator (FT of GF)
op=bsxfun(@times, omega, op); % Magic operator

E=real(-fftTimes(op, rho));
B=real(fftCross(op, J));

figure(1); axis equal;
hold on;
quiver3(xx,yy,zz,E(:,:,:,1),E(:,:,:,2),E(:,:,:,3));
quiver3(xx,yy,zz,B(:,:,:,1),B(:,:,:,2),B(:,:,:,3));
end