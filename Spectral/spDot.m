function w=spDot(filter, u)
% Returns the dot product in the Fourier space.
u_hat(:,:,:,1)=fftn(u(:,:,:,1));
u_hat(:,:,:,2)=fftn(u(:,:,:,2));
u_hat(:,:,:,3)=fftn(u(:,:,:,3));
w_hat=dot(filter, u_hat, 4);
w=ifftn(w_hat);
end