function [] = blocki(n, m)
adjx=[2 3; 3 4; 5 6; 6 7];
adjy=[3 1; 1 6];

[topo,net,RL,TB]=ddtopo(adjx,adjy);
pos=ddpatches(topo);

figure(1); clf; hold on;

x=linspace(-1,1,n);
y=linspace(-1,1,m);
[xx,yy]=ndgrid(x,y);
zq=xx+1i*yy;

for i=1:length(pos)
    zz=zq+pos(i);
    surf(real(zz), imag(zz), 0*real(zz));
    colormap([1.0,0.4,0.0]);
end
hold off;
end