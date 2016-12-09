function [img,sub,bis,sup] = img_sq(n,flag)
%IMG_SQ  Partition a square grid into sub-domains above, below and on the
%        bisector, identifying their respective image-points.
%
% [img,sub,bis,sup] = img_sq(n)
%     Identify bisector with main diagonal [notionally y=x].
%
% [img,sub,bis,sup] = img_sq(n,1)
%     Identify bisector with reverse diagonal [notionally y = -x].
%
% img - indices of image-points [integer, length n^2, reordering of 1:n^2]
% sub -   below-bisector points [logical, length n^2]
% bis -         bisector points [logical, length n^2]
% sup -   above-bisector points [logical, length n^2]

    if nargin < 2,  flag = 0;  end

    ind = (1:n^2)';
    Ind = reshape(ind,n,n);
    
    if flag,  Ind = fliplr(Ind);  end
    
    Img = Ind';
    Sub = triu(Ind,1);
    Bis = diagmat(Ind);
    
    if flag  % perform reverse-flip
        Img = fliplr(Img);
        Sub = fliplr(Sub);
        Bis = fliplr(Bis);
    end
    
    img = Img(:);
    sub = logical(Sub(:));
    bis = logical(Bis(:));
    sup = ~(sub | bis);

end
