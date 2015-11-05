function []=pdesolver(n)
[D,x]=chebD(n+2);
D2=D*D; D2=D2(2:end-1,2:end-1);

[xx, yy]=meshgrid(x);
F=cos(pi*xx).*cos(pi*yy);
F=F(2:end-1,2:end-1);

I=eye(n);
X=zeros(n);
[Q, U]=schur(D2);
FF=Q'*F*Q;
for i=n:-1:1
    X(:, i)=(U+U(i,i)*I)\(FF(:,i)-X(:,i+1:end)*U(i,i+1:end)');
end
X=Q*X*Q';

disp(norm(D2*X+X*D2'-F, 'fro'))
uu=zeros(n+2);
uu(2:end-1,2:end-1)=X;
figure(2);
surf(xx,yy,uu);
end