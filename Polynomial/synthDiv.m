function [Q,R] = synthDiv(P, a)
n=length(P);
Q(n-1)=P(n);
for i=n-2:-1:1
    Q(i)=Q(i+1)*a+P(i+1);
end
R=Q(1)*a+P(1);
end
