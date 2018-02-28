% Sample gird
N=32;
x=linspace(-1,1,N)';
y=linspace(-1,1,N)';
[xx,yy]=ndgrid(x,y);



% Vertices (EN, WN, ES, WS)
Z=[(1+1i), (-1+1i)*4, (1-1i), (-2-1i), (1+4i), (-3+6i)];
% Radius of curvautre (E, W, N, S)
R=[inf,23,7,inf; inf,inf,inf,7];

quad=[1,2,3,4; 5,6,1,2];


figure(1); clf;
for k=1:size(quad,1)
    z=Z(quad(k,:));
    r=R(k,:);
    
    F = curvedquad(z,r);
    [jac,G11,G12,G22] = diffgeom(F,x,y); 
    zz = F(xx,yy);
    
    subplot(2,2,1); 
    surf(real(zz), imag(zz), jac); hold on;
    title('Jacobian');
    colormap(jet(256));
    %shading interp;
    axis square;
    %camlight;
    view(2);

    subplot(2,2,2);
    surf(real(zz), imag(zz), G11); hold on;
    title('Metric G11');
    colormap(jet(256));
    %shading interp;
    axis square;
    %camlight;
    view(2);

    subplot(2,2,3);
    surf(real(zz), imag(zz), G12); hold on;
    title('Metric G12');
    colormap(jet(256));
    %shading interp;
    axis square;
    %camlight;
    view(2);

    subplot(2,2,4);
    surf(real(zz), imag(zz), G22); hold on;
    title('Metric G22');
    colormap(jet(256));
    %shading interp;
    axis square;
    %camlight;
    view(2);
end