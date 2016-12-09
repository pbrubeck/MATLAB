function r = spdensity(A)
%SPDENSITY  Proportion of non-zero entries in a sparse matrix.

    if issparse(A),
        r = length(find(A))/length(A(:));
    else
        r = 1;
    end

end

