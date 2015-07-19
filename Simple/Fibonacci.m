function fib = Fibonacci(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fib=zeros(n,1);
fib(2)=1;
for i=3:n
    fib(i)=fib(i-1)+fib(i-2);
end

end

