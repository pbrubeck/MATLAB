function [J] = barinterp(x1,x2)
N1=length(x1);
N2=length(x2);
lam=repmat(x1(:),1,N1);
lam=lam-lam'+eye(N1);
lam=1./prod(lam,1);
J=repmat(x2(:),1,N1)-repmat(x1(:)',N2,1);
[ii,jj]=find(abs(J)<10*eps);
J=repmat(lam,N2,1)./J;
J(ii,:)=0; J(ii+N2*(jj-1))=1;
J=J.*repmat(1./sum(J,2),1,N1);
end