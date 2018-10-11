function [u] = hank(B,jkrm,u)
n=size(u,2);
m=fftshift(-n/2:n/2-1);
u=fft(B*u,[],2);
for j=1:size(jkrm,3)
    u(:,j)=jkrm(:,:,j)*u(:,j);
end
u=ifft(u*diag(1i.^m),[],2);
end