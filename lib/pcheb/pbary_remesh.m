function [C,x2] = pbary_remesh(N1,N2,xel1,xel2)
%PBARY_REMESH  Interpolate from one PCHEB grid to another.
% 
% [C,x2] = pbary_remesh(N1,N2,xel)
%    Interpolate from X1 = PBARY(N1,xel) to X2 = PBARY(N2,xel).
%    C is of size n2-by-n1, where n1=length(x1) and n2=length(x2).
% 
% [C,x2] = pbary_remesh(N1,N2,xel1,xel2)
%     Interpolate from X1 = PBARY(N1,xel1) to X2 = PBARY(N2,xel2).

    if nargin < 4,  xel2 = xel1;  end

    [~,x2,C] = pcheb(N1,xel1, pcheb(N2,xel2,[]) );

end
