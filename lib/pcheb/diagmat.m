function A = diagmat(A)
%DIAGMAT Square matrix with off-diagonal elements zeroed.

    A = diag(diag(A));

end
