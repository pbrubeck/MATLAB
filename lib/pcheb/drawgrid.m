function drawgrid(x,y,flag)
%DRAWGRID  Show domain decomposition (x,y), where vectors X and Y may be
%          either collocation points or element boundaries.

    if nargin < 3,  flag = 0;  end
    
    switch flag
        case 1
            drawgrid_sq(x)
        case 2
            drawgrid_rect(x,y),  axis equal
        case 3
            drawgrid_lwise(x,y)
        otherwise
            drawgrid_rect(x,y)
    end
    
end
