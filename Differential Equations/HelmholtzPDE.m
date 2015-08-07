function A=HelmholtzPDE(n)
% Solves the Helmholtz problem in two dimensions using Chebyshev methods
[Del2, xx, yy]=chebLaplacian(n);
k=9;
op=Del2+k^2*eye((n-2)^2);
[L, U]=lu(op);

t=0; dt=0.01;
nframes=1000;
A=zeros(n-2,n-2);
figure(1)
h=surf(xx, yy, A, 'EdgeColor', 'none');
colormap(jet);
for i=0:nframes
f=exp(-30*((xx-0.5*cos(t)).^2+(yy-0.5*sin(t)).^2));
A=reshape(U\(L\f(:)), size(f));
set(h, 'ZData', A);
drawnow;
t=t+dt;
end

end

