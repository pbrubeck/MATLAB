function [] = waveletDenoise()
X=imread('Lenna.png');
X=im2double(X);

% Gaussian noise
sigma=0.08;
Y=imnoise(X,'gaussian',0,sigma^2);

% Haar wavelet transform
N=wmaxlev(min(size(X(:,:,1))),'haar');
[C,S]=wavedec2(Y,N,'haar');
% Kill lower coefficients
CT=C; CT(abs(C)<3*sigma)=0; 
% Keep coarse coefficients
k=3*prod(S(1,:));
CT(1:k)=C(1:k);
% Inverse Haar wavelet transform
Z=waverec2(CT,S,'haar');
p=psnr(Z,X);

% Show images
figure(4); 
imagesc([X,Y,Z]); colormap(gray); 
title(sprintf('Wavelets PSNR = %f', p));
end