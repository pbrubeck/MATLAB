function [] = imedges()
A=im2double(rgb2gray(imread('lena.tiff')));
C=fspecial('sobel');
G=abs(conv2(A, C+1i*C', 'same'));
G=im2bw(G/max(G(:)), 0.13);
imagesc(G); colormap(gray); axis equal off;
end