function [AA,BB,Q,Z,order,start,csize,lambda] = clustered_qz(A,B,epscluster)

%CLUSTERED_QZ  QZ with clustering 
%
% [AA,BB,Q,Z,order,start,csize,lambda] = CLUSTERED_QZ(A,B,epscluster)

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

[tmpS0,tmpS1,tmpQ,tmpZ] = qz(A,B);
tmplambda = diag(tmpS1)./diag(tmpS0);
[order,start,csize,lambda] = clusters(tmplambda,epscluster);
n = length(order);
reverse = order(n:-1:1);
invord = zeros(1,n);
for j = 1:length(reverse)
    invord(reverse(j)) = j;
end
[AA,BB,Q,Z] = ordqz(tmpS0,tmpS1,tmpQ,tmpZ,invord);

