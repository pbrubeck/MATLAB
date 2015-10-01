function [E, B] = Maxwell(rho, J)
% Solves Maxwell's (static) equations to find E and B
a=-pi; b=pi;
n=size(rho,1);
N=[n,n,n];
P=b-a;
th=2i*pi/P;
gv=a+P/n*(0:n-1);
[xx,yy,zz]=meshgrid(gv);

i=th*[0:n/2, -n/2+1:-1];
j=th*[0:n/2, -n/2+1:-1];
k=th*[0:n/2, -n/2+1:-1];
[ii,jj,kk]=meshgrid(i, j, k);
omega=cat(4, ii, jj, kk); % Frequency space
[i2,j2,k2]=meshgrid(i.^2, j.^2, k.^2);
D2=i2+j2+k2; % Laplacian operator

px=exp(-a*i); py=exp(-a*j); pz=exp(-a*k);
phase=reshape(kron(kron(px,py),pz), N(1:3));
op=phase./D2; op(1,1,1)=0; % Newtonian potential operator
op=bsxfun(@times, omega, op); % Magic operator


E=real(-spTimes(op, rho));
B=real(spCross(op, J));

figure(1); axis equal;
hold on;
quiver3(xx,yy,zz,E(:,:,:,1),E(:,:,:,2),E(:,:,:,3));
quiver3(xx,yy,zz,B(:,:,:,1),B(:,:,:,2),B(:,:,:,3));
end