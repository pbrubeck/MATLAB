function out1=Example1(v)

n = size(v,1); J=2:n-1;
out1=[2*v(1,:)-v(2,:); 2*v(J,:)-v(J-1,:)-v(J+1,:); -v(n-1,:)+2*v(n,:)];
out1=n^2*out1;


return
