function [] = invfilter()
% Load image
im0=rgb2gray(imread('lena.tiff'));
im0=double(im0)/max(double(im0(:)));

filter=fspecial('laplacian');

% Apply Fourier filter -> Convolution
filter_hat=fft2(filter, size(im0,1), size(im0,2));
filter_hat=filter_hat+(filter_hat==0); % non-singular filter
im1=ifft2(fft2(im0).*filter_hat);

figure(1);
imagesc(im1); colormap(gray); colorbar(); axis equal;

% Invert filter -> Deconvolution
imhat=fft2(im1)./filter_hat; 
im2=ifft2(imhat);

figure(2);
imagesc(im2); colormap(gray); colorbar(); axis equal;
end

