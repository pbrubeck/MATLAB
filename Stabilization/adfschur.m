function [outputArg1,outputArg2] = adfschur(n,nel,nu)

[D,x,w]=legD(n);



nex=nel;
ntot=n*nel;
CFL=0.001;
%CFL=inf;

vx=1*ones(1,nel);


rd=[1];

mask=ones(n,nel,1);
mask(rd)=0;


cid=repmat((1:2)',1,nel)+repmat((0:nel-1)*(2-1),2,1);


[D,zhat]=legD(n);
Vhat=VandermondeLeg(zhat);
Bhat=inv(Vhat*Vhat');
Ahat=D'*Bhat*D;
Chat=Bhat*D;
HPF=Bhat*(bubfilt(zhat)-eye(n));

xc=linspace(-1,1,nex+1);
h=diff(xc)/2;

x=zhat*h+repmat(conv(xc,[1,1]/2,'valid'),n,1);
dx=min(diff(zhat))*min(h(:));
dt=CFL*dx;

aa=zeros(n,n,nel);
ia=zeros(n,n,nel);
ja=zeros(n,n,nel);
for e=1:nel
    aa(:,:,e)=(nu./h(e))*Ahat+vx(e)*Chat+(h(e)/dt)*HPF;
    [ia(:,:,e),ja(:,:,e)]=ndgrid(1+(e-1)*(n-1):1+e*(n-1));
end
A=sparse(ia(:),ja(:),aa(:));

rd=[1,size(A,1)];
A(rd,:)=0;
A(:,rd)=0;
A(rd,rd)=speye(length(rd));

j=n:n-1:1+(nel-1)*(n-1);
i=setdiff(1:size(A,1),j);

S=A(j,j)-A(j,i)*(A(i,i)\A(i,j));
imagesc(full(S));
end

