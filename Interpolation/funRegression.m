function c = funRegression(x, y, F)
% Returns the coeficients of the linear combination of the functions on F
% that best fit the data.
A=transpose(F(x));
c=(A'*A)\(A'*y(:));
end