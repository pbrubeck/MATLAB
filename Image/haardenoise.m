function [Xhaar, C, CT, S, lvl] = haardenoise(X, sigma)
% Perform Haar Wavelet Transform
lvl=wmaxlev(min(size(X(:,:,1))),'haar');
[C,S]=wavedec2(X,lvl,'haar');
% Kill the low coefficients in last three detail levels
CT=C;
CT((1:end>prod(S(end-3,:))) & (abs(C)<1*sigma))=0;
CT((1:end>prod(S(end-2,:))) & (abs(C)<2*sigma))=0;
CT((1:end>prod(S(end-1,:))) & (abs(C)<3*sigma))=0;
% Inverse Haar Wavelet Transform
Xhaar=waverec2(CT,S,'haar');
end