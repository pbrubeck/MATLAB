function [SA] = schwarz1d(n,no,A,lbc,rbc)
ns=n+2*no;
SA=zeros(ns,ns,size(A,3));

i1=1:no+1;
i0=i1(end):i1(end)+n-1;
i2=i0(end):i0(end)+no;

j1=1:no+1;
j0=j1(end)+1:j1(end)+n;
j2=j0(end)+1:j0(end)+1+no;

SA(i0,i0,:)=A(j0,j0,:);
for k=1:size(A,3)
    if(lbc==0) % Neumann      
        SA(i1(1:no),i1(1:no),k)=eye(no);
    elseif(lbc==1) % Dirichlet
        SA(i1(end),:,k)=0;
        SA(:,i1(end),k)=0;
        SA(i1,i1,k)=eye(no+1);
    elseif(lbc==2) % Overlap
        SA(i1,i1,k)=SA(i1,i1,k)+A(j1,j1,k);
    end    
    if(rbc==0) % Neumann
        SA(i2(end-no+1:end),i2(end-no+1:end),k)=eye(no);
    elseif(rbc==1) % Dirichlet
        SA(i2(1),:,k)=0;
        SA(:,i2(1),k)=0;
        SA(i2,i2,k)=eye(no+1);
    elseif(rbc==2) % Overlap
        SA(i2,i2,k)=SA(i2,i2,k)+A(j2,j2,k);
    end
    
end
end