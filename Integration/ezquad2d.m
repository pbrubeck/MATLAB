function [xx,yy,ww] = ezquad2d(x,y,wx,wy)
% use I=ww*f(xx,yy);
[xx,yy]=meshgrid(x,y);
ww=kron(wx,wy);
xx=xx(:);
yy=yy(:);
ww=ww(:).';
end