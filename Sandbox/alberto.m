function [] = alberto( n )
% input data
W=1+rand(n,n);
f1=1+rand(n,1);
f2=1+rand(n,1);

% indices
m=2*n+2;
j1=1:n;
j2=n+1:2*n;
j11=1:m:n*m; % first half diagonal
j22=1+n*m:m:2*n*m; % second half diagonal


% initial values
x=1+rand(n,1);
y=1+rand(n,1);
z=rand(1,1);
u=[x;y;z];
J=zeros(2*n+1);
J(:,end)=[-f1;-f2;-1];
err=1;
i=1;
while(err>1E-10)
    x=u(1:n);
    y=u(n+1:2*n);
    z=u(2*n+1);
    
    Wx=W*x;
    Wy=W*y;
    
    J(j11)=Wy;
    J(j22)=Wx;
    J(j1,j2)=bsxfun(@times, x, W);
    J(j2,j1)=bsxfun(@times, y, W);
    J(end,j1)=Wy;
    J(end,j2)=x'*W;
    
    f=[x.*Wy-z*f1; y.*Wx-z*f2; x'*Wy-z];
    du=J\f;
    u=u-du;
    err=norm(du,'inf');
    i=i+1;
end
display(i);
end

