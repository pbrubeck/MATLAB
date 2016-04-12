function [a] = TaylorSeries(f, x0, n)
x=linspace(-1,1,2*n);
y=f(x-x0);
a=polyRegression(x,y,n);
end

