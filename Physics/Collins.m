function [G] = Collins(E, xx, yy, ww)
% Parallel Gauss Transform
G=gpuArray(zeros(size(E)));
xx=gpuArray(xx);
yy=gpuArray(yy);
ww=gpuArray(ww);
k=1;
A=[1,0;0,1];
B=[1,0;0,1];
C=[1,0;0,1];
D=[1,0;0,1];

r1=[xx, yy];
T=dot(r1*(A\B), r1, 2);

parfor n=1:numel(E)
    r2=[xx(n);yy(n)];
    K=r2
    G(n)=ww*(E.*exp(-1i*k/2*(T)));
end
end