function [] = imedges()
%EDGES Summary of this function goes here
%   Detailed explanation goes here
im=rgb2gray(imread('Lenna.png'));
C=fspecial('sobel');

C1_hat=fft2(C, size(im,1), size(im,2));
C2_hat=fft2(C', size(im,1), size(im,2));
im_hat=fft2(im);

G1=ifft2(im_hat.*C1_hat);
G2=ifft2(im_hat.*C2_hat);

G=sqrt(G1.^2+G2.^2);
imagesc(G); colormap(gray); colorbar(); axis equal;
end

