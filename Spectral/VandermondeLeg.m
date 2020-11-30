function [ V ] = VandermondeLeg(x,m)
if(nargin<2), m=length(x); end
V=zeros(length(x),m);
for j=1:m
    V(:,j) = LegendreP([zeros(j-1,1);1], x);
    V(:,j) = sqrt((2*j-1)/2) * V(:,j);
end
end