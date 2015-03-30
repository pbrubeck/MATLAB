function poly = polyRegression(data, deg)
A=ones(length(data),deg+1);
for i=1:length(data)
    for j=1:deg
        A(i,j+1)=A(i,j)*(i-1);
    end
end
poly=(A'*A)\(A'*data);
end