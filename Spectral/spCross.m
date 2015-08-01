function w=spCross(filter, u)
% Returns the cross product in the Fourier space.
u_hat(:,:,:,1)=fftn(u(:,:,:,1));
u_hat(:,:,:,2)=fftn(u(:,:,:,2));
u_hat(:,:,:,3)=fftn(u(:,:,:,3));
w_hat=cross(filter, u_hat);
w(:,:,:,1)=ifftn(w_hat(:,:,:,1));
w(:,:,:,2)=ifftn(w_hat(:,:,:,2));
w(:,:,:,3)=ifftn(w_hat(:,:,:,3));
end