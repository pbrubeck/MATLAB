function out1=Example2(v,flag)

n = 200;
J=2:n-1;

if nargin <2

   out1=[2*v(1,:)-v(2,:); 2*v(J,:)-v(J-1,:)-v(J+1,:); -v(n-1,:)+2*v(n,:)];

else

   switch(flag)

   case  'dimension'
      out1 = n;      %%% dimension of the matrix operation

   case 'preconditioner'
      e = ones(n,1); 
      L=spdiags([-e,e],-1:0,n,n);
      U=L';

      out1 = U\(L\v);

   case 'L'
      e = ones(n,1); 
      L=spdiags([-e,e],-1:0,n,n);

      out1 = L\v;

   case 'U'
      e = ones(n,1); 
      U=spdiags([e,-e],0:1,n,n);

      out1 = U\v;


   otherwise
      error(['Unknown flag ''' flag '''.']);
   end
end

return
