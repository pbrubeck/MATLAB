function quadplot(x,y,U,xx,yy,Cx,Cy,flag)
%QUADPLOT  Draw contour plot, two mesh plots, and a grid-plot.

    if nargin < 8,  flag = 0;  end

    UU = expand2d(U,Cx,Cy);

    subplot(2,2,1),  mycontour(xx,yy,UU,flag)

    subplot(2,2,3),  drawgrid(x,y,flag)

    subplot(2,2,2),  mesh(xx,yy,UU),  xlabel x,  ylabel y
    
    xlim([xx(1) xx(end)])
    ylim([yy(1) yy(end)])

    subplot(2,2,4),  mesh(x,y,U),  xlabel x,  ylabel y
    
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])

end

