function drawgrid_rect(x,y)
%DRAWGRID_RECT  Show domain decomposition for SOLVE_LIN_RECT.

    x0 = x(1);  x1 = x(end);
    y0 = y(1);  y1 = y(end);
    
    for i = 1:length(x),
        xi = x(i);  plot([xi,xi], [y0,y1], 'k-'),  hold on
    end
    
    for i = 1:length(y),
        yi = y(i);  plot([x0,x1], [yi,yi], 'k-'),  hold on
    end
    
    hold off

end
