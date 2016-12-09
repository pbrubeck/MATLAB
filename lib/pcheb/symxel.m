function xel = symxel(L,rho)
%SYMXEL  Symmetrised element boundaries on [-L,L], including zero.
%        Vector 0 < RHO < 1 specifies non-terminal elements in (-L,L).

    r = rho(:)';  r = sort(r(r>0 & r<1));
    
    xel = abs(L) * [-1,-fliplr(r),0,r,1];

end

