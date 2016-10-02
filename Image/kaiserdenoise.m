function [Xkaiser, betopt, pmax] = kaiserdenoise(X, X0)
% Finds the Kaiser window that maximizes psnr
beta=0:48;
pmax=-inf; p=zeros(size(beta));
[f1,f2]=freqspace(9,'meshgrid');
Hd=hypot(f1, f2)<0.25;
for i=1:length(beta);
    h=fwind1(Hd,kaiser(9,beta(i)));
    Z=imfilter(X,h,'replicate');
    Z=Z/max(Z(:));
    p(i)=psnr(Z, X0);
    if(p(i)>pmax)
        pmax=p(i);
        betopt=beta(i);
        Xkaiser=Z;
    end
end
end

