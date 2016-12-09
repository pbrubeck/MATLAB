function drawgrid_lwise(x,y)
%DRAWGRID_LWISE  Show domain decomposition for SOLVE_LIN_LWISE codes.

    L = max(x);
    
    e = 0.05;  % gap 
        
    for i = 1:length(x),
        xi = x(i);
        ymin = -L*(xi <= 0);
        plot([xi,xi], [ymin,L], 'k-'),  hold on
    end
    
    for i = 1:length(y),
        yi = y(i);
        xmax = L*(yi >= 0);
        plot([-L,xmax], [yi,yi], 'k-'),  hold on
    end
    
    hold off,  axis((1+e)*L*[-1 1 -1 1])
    
    axis_sq

end
