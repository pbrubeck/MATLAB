function [ A ] = complex2rgb(U, p)
if(nargin==1)
    p=0;
end
% Square modulus defines value and saturation
I=abs(U).^2;
V=p+(1-p)*mat2gray(I);
V=cat(3, V, V, V);
% Phase defines hue
hue=angle(U)/(2*pi);
hue=hue-floor(hue);
map=hsv(256);
H=ind2rgb(uint8(255*hue), map);
A=H.*V;
end

