
N=8;
Ex=[20,40];
Ey=[20,40];
nu=[-25,-200,-800];
%nu=-10.^(0:2:4);
ifcrs=1;
ifdeal=1; 

opts=ones(8,5);
opts(1,1:2)=[-1,0];
opts(2,1:2)=[-1,0];
opts(3,1:2)=[0,0];
opts(4,1:2)=[0,0];
opts(5,1:2)=[1,0];
opts(6,1:2)=[0,0];
opts(7,1:2)=[1,1];
opts(8,1:2)=[1,1];

opts(:,3)=ifcrs;
opts(:,4)=repmat([0;1],size(opts,1)/2,1);
opts(:,5)=ifdeal;

[nu,Ey,Ex]=ndgrid(nu,Ey,Ex);
disp([Ex(:), Ey(:), nu(:)]);

iter=zeros(length(nu),size(opts,1));
for j=4:size(opts,1)
    for i=1:numel(nu)
        [relres,it,resvec]=adf2(N,[Ex(i),Ey(i)],nu(i),opts(j,:));
        iter(i,j)=it;
    end
    disp(iter(:,1:j));
end
iter=reshape(iter,[],size(opts,1));
display(iter);