function drawgrid_lrefine(xel,rho)
%DRAWGRID_LREFINE  Show domain decomposition for SOLVE_LIN_SQ_LREFINE.

    x0 = xel(1);  x1 = xel(2);  xelc = x0 + (x1-x0)*[ 0; rho(:); 1 ];

    drawgrid_sq(xel),   hold on
    drawgrid_sq(xelc),  hold off

end

