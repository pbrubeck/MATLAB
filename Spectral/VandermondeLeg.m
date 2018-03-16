function [ V ] = VandermondeLeg(x)
V=zeros(length(x));
for j=1:length(x)
    V(:,j)=sqrt((2*j-1)/2)*LegendreP([zeros(j-1,1);1], x);
end
end