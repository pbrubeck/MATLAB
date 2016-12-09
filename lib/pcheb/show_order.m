function show_order(dof,Nx,Ny)
%SHOW_ORDER  Pretty-print polynomial order of spectral-element scheme.

    if nargin < 3,  Ny = Nx;  end
    
    str1 = ['dof = ',int2str(dof)];
    
    if length(Nx)==length(Ny) && all(Nx==Ny),
        str2 = ['N = ',int2str(Nx)];
    else
        str2 = ['Nx = ',int2str(Nx),', Ny = ',int2str(Ny)];
    end
        
    disp([str1,'  (',str2,')'])

end

