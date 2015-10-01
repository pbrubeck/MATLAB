function w=spTimes(filter, u)
% Returns element wise multiplication in the Fourier space.
u_hat=fftn(u);
w_hat=bsxfun(@times, filter, u_hat);
w(:,:,:,1)=ifftn(w_hat(:,:,:,1));
w(:,:,:,2)=ifftn(w_hat(:,:,:,2));
w(:,:,:,3)=ifftn(w_hat(:,:,:,3));
end