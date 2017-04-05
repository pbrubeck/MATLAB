function [] = riemmanhyp( N )
[D,x]=lagD(N);
[V,L]=eig(-1i*D,'vector');
W1=V*diag(1./(1-exp(-1i*L)))/V;
W2=V*diag(1-exp(-1i*L))/V;
H=-1i*W1*(diag(x)*D+D*diag(x))*W2;

kd=2:N;
Psi=zeros(N);
E=zeros(N,1);
[Psi(kd,kd),E(kd)]=eig(H(kd,kd),'vector');
[~,id]=sort(abs(imag(E)));

Psi=Psi(:,id);
E=E(id);

disp(E);
plot(x,Psi);
end