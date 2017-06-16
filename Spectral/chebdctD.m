function v=chebdctD(u)
n=length(u);
N=2*n-2;
uhat=ifft(u,N,'symmetric');
uhat=uhat(1:n);

A=triu(repmat(0:2:2*n-2,n,1));
A(1:2:end,1:2:end)=0;
A(2:2:end,2:2:end)=0;
uhat=A*uhat;

v=N*ifft(uhat,N,'symmetric');
v=v(1:n);
end