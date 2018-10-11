function [jkrm] = hankbess(n,r)
% n : angular modes
% r : quadrature nodes [0,L]
kr=r(:)*r(:)';
jkrm=zeros(length(r),length(r),n);
m=fftshift(-n/2:n/2-1);
for i=1:n
jkrm(:,:,i)=besselj(m(i),kr);
end
end