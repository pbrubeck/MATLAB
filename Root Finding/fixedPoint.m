function xi = fixedPoint(f, xi)
%Solves x=f(x) for x
ii=0;
prev=0;
while(ii==0 || abs(xi-prev)>(1E-15)*xi)
    fprintf('i=%d \t xi=%f \t err=%f \n', ii, xi, 100*abs((xi-prev)/xi));
    prev=xi;
    xi=f(xi);
    ii=ii+1;
end
end

