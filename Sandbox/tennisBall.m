function [] = tennisBall(N)
figure(1);
[h,ph,th]=sphericalPlot(2*N,N);

% A=0.5;
% gamma=0.35;
% t=cos(2*ph);
% u=sign(t).*(A*abs(t).^gamma);
% C=sign(2/pi*th-1-u);

C1=sign((sin(th).*cos(ph)).^2-0.5);
C2=-sign((sin(th).*sin(ph)).^2-0.5);
idx=th<pi/2;
C1(idx)=C2(idx);

alpha(0.6);
set(h,'CData',C1);
shading interp;
colormap([1 0 0; 0 0 1])
end