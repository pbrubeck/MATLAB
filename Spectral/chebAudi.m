function [D, x] = chebAudi(N, m)
x=cos(pi*(0:N-1)/(N-1))';

z = ainit(x,m+2); % x0 = 1
y = sin((N-1)*acos(z)); % y = f(x)
y.k = y.k - 1; % reduce order of numerator
D=zeros(N);
for i=1:N
    q = y./(adiff(y,1).*(z-x(i))); % result (ignore entries of q.c beyond "nan")
    D(:,i)=q{m}/(1+(i==N)+(i==1));
end
end