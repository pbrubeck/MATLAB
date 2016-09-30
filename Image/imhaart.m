function [] = imhaart()
% Black and white image
X=imread('Lenna.png');
X=double(X)/255;

% Gaussian noise
sigma=0.08;
w=sigma*randn(size(X));
Y=max(min(X+w,1),0);

% Haar wavelet transform
N=wmaxlev(min(size(X(:,:,1))),'haar');
[C,S]=wavedec2(Y,N,'haar');

% Kill lower coefficients
tol=3*sigma;
CT=C;
CT(abs(C)<tol)=0;
CT(1:3*S(1,1)*S(1,2))=C(1:3*S(1,1)*S(1,2));

% Inverse Haar wavelet transform
Z=waverec2(CT,S,'haar');

% Show images
figure(1); 
imagesc([X,Y,Z]);
axis equal off;
end