function [Xkaiser, Xhaar]=imdenoise(filename, sigma)
% Attemps to denoise an image via Kaiser window and Haar wavelet transform.
% Examples: 
%   [Xkaiser, Xhaar]=imdenoise('Lenna.png', 0.08);
%   [Xkaiser, Xhaar]=imdenoise('airplane.png', 0.1);

% The best Kaiser window is found with previous knowledge of the original
% image in order to compare the accuracy of the wavelet filtering.

%% Add gaussian noise to image
close all;
X0=im2double(imread(filename));
X=imnoise(X0,'gaussian',0,sigma^2);

% Show the noisy image
figure(1); imshowpair(X0, X, 'montage');
title(sprintf('Original | Gaussian noise psnr = %.4f', psnr(X,X0)));

%% Find the best Kaiser window
[Xkaiser, betopt, pmax]=kaiserdenoise(X,X0);

% Show reconstructed image
figure(2); imshowpair(X0, Xkaiser , 'montage');
title(sprintf('Original | Kaiser 9x9 window beta = %.1f psnr = %.4f', betopt, pmax));

%% Haar wavelet transform denoising
[Xhaar, C, CT, S, lvl]=haardenoise(X,sigma);
% Gaussian blur
Xhaar=imfilter(Xhaar, fspecial('gaussian',3));

% Show reconstructed image
figure(3); imshowpair(X0, Xhaar, 'montage');
title(sprintf('Original | Haar Wavelet Transform psnr = %.4f', psnr(Xhaar,X0)));

% Plot wavelet coefficients
figure(4); plotwavelet2(C(1:size(X,3):end),S(:,1:2),lvl,'haar',255,'square'); colormap(gray); title('Noisy wavelet coefficients (R channel)');
figure(5); plotwavelet2(CT(1:size(X,3):end),S(:,1:2),lvl,'haar',255,'square'); colormap(gray); title('Filtered wavelet coefficients (R channel)');