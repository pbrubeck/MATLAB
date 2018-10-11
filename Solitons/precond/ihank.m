function [u] = ihank(B,jkrm,u)
n=size(u,2);
m=fftshift(-n/2:n/2-1);
u=ifft(B*u,[],2)*diag(1i.^-m);
for j=1:size(jkrm,3)
    u(:,j)=jkrm(:,:,j)*u(:,j);
end
u=fft(u,[],2);
end