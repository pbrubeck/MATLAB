function out=digitalComb(x)
n=6;
A=de2bi((0:2^n-1)', n);
b=[106,80; 193,120; 303,280; 50,40; 140,100; 243,200];
z=A*b;
out=zeros([numel(x), n]);
for i=1:numel(x)
    w=z(:,2);
    w(z(:,1)<x(i))=inf;
    [~, idx]=min(w);
    out(i,:)=A(idx(1), :);
end
end