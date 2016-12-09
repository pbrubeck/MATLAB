function drawgrid_refine(xel,rho)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



    x0 = xel(1);  x1 = xel(2);  xend = xel(end);  dx = x1 - x0;

    rnge = [x0,xend];

    for i = 1:length(xel),
       xi = xel(i);
       plot([xi,xi], rnge, 'k-'),  hold on
       plot(rnge, [xi,xi], 'k-'),  hold on
    end

    rnge = [x0,x1];

    for i = 1:length(rho),
       xi = x0 + dx*rho(i);
       plot([xi,xi], rnge, 'k-'),  hold on
       plot(rnge, [xi,xi], 'k-'),  hold on
    end


    hold off
    
    axis_sq


end

