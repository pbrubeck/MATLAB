function w=spTimes(op, u)
% Returns element-wise multiplication in the Fourier space.
if ndims(u)==4
    u_hat(:,:,:,1)=fftn(u(:,:,:,1));
    u_hat(:,:,:,2)=fftn(u(:,:,:,2));
    u_hat(:,:,:,3)=fftn(u(:,:,:,3));
else
    u_hat=fftn(u);
end
w_hat=bsxfun(@times, op, u_hat);
w(:,:,:,1)=ifftn(w_hat(:,:,:,1));
w(:,:,:,2)=ifftn(w_hat(:,:,:,2));
w(:,:,:,3)=ifftn(w_hat(:,:,:,3));
end

