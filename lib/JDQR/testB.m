function test

clear all
close all hidden
rand('seed',0);


method='jdqr';

nselect=10;
SHOW=2;
sigma=0.0;
tol=1.e-12;

J=[0,1]; 

SelectExample(10)=1;

SelectExample=[SelectExample,zeros(1,20)];

for k=1:length(J), TS=J(k);


  for Ex=1:length(SelectExample)
  rand('seed',0)
  if SelectExample(Ex)


    fprintf('\n\n=========================')
    fprintf('===================================')
    fprintf('\nExample %i',Ex)
    fprintf('\nSolution method:    %s\n\n',method)
    if Ex>7, SHOW0=SHOW; else, SHOW0=max(0,SHOW-1); end

    A=Example(Ex);
    options=struct('Tol',tol,'Disp',SHOW0,'TestSpace',TS,'Schur','no'); 
    [Xeig,Lambda]=feval(method,A,nselect,sigma,options); 
 
    nrm=norm(A*Xeig-Xeig*Lambda);
    if nrm> tol
       fprintf('\nnorm(%s*X-X*Lambda): %0.5g\n','A',nrm), pause
    end

    if Ex>7, drawnow, figure, end

    end
  end
end

close

return


%=========================================================
function A=Example(Ex)

  jj=sqrt(-1);

  switch Ex
  case 1
    A=zeros(5,5); A(:,5)=rand(5,1);     %% 1
  case 2
    A=zeros(5,5);                       %% 2
  case 3
    A=zeros(5,5); A(1:3,4:5)=rand(3,2); %% 3
  case 4
    A=rand(3,3);                        %% 4
  case 5
    A=rand(5,5)+jj*rand(5,5);           %% 5
  case 6
    A=rand(8,8)+jj*rand(8,8);
    A=triu(A,1); A(1,1)=1; A(2,2)=1; 
    A(4,4)=rand+jj*rand;                %% 6
  case 7
    n=17;
    j= 1138264952; 
    rand('seed',j)
    A=rand(n,n); A=triu(A);             %% 7
  case 8
    n=100; e=ones(n,1);
    A=spdiags([e -2*e e],-1:1,n,n);     %% 8
  case 9
    n=100; e=ones(n,1);
    A=spdiags([e -2*e 1.2*e],-1:1,n,n);  %% 9
  case 10
    rand('seed',164446964) 
    n=100; e=ones(n,1); o=zeros(n,1);
    A=spdiags([e o -2*e e],-2:1,n,n);   %% 10
  end

return
