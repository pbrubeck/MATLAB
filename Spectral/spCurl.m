function curl = spCurl(u)
% Calculates the spectral curl of a vector field, actual curl
% must be scaled by 2pi/L.
N=size(u);
i=[0:N(1)/2-1, 0, -N(1)/2+1:-1];
j=[0:N(2)/2-1, 0, -N(2)/2+1:-1];
k=[0:N(3)/2-1, 0, -N(3)/2+1:-1];
[ii,jj,kk]=meshgrid(i, j, k);
omega=cat(4, ii, jj, kk);
u_hat(:,:,:,1)=fftn(u(:,:,:,1));
u_hat(:,:,:,2)=fftn(u(:,:,:,2));
u_hat(:,:,:,3)=fftn(u(:,:,:,3));
curl_hat=1i*cross(omega, u_hat);
curl(:,:,:,1)=ifftn(curl_hat(:,:,:,1));
curl(:,:,:,2)=ifftn(curl_hat(:,:,:,2));
curl(:,:,:,3)=ifftn(curl_hat(:,:,:,3));
end