function mycontour(x,y,U,flag)
%MYCONTOUR  Produce a customised contour plot.

    if nargin < 4,  flag = 0;  end

    contour(x,y,U), xlabel x,  ylabel y
    
    switch flag
        case 1  % square mode
            axis_sq
        case 2  % equal mode
            axis equal
            axis([x(1) x(end) y(1) y(end)])
        case 3  % L-wise mode
            axis_sq
            L = max(max(x));
            hold on
            plot([0,L],[ 0,0],'k')  % draw wall at y==0 (horizontal)
            plot([0,0],[-L,0],'k')  % draw wall at x==0 (vertical)        
            hold off
    end
    
end
