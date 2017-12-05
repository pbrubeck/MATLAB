function [topo,net,RL,TB] = ddtopo( adjx, adjy )
% Domain Decomposition Topology
% 2D Tensor product domains

% adjx : Horizontal interfaces [east  domain, west  domain]
% adjy : Vertical interfaces   [north domain, south domain]

adj=[adjx; adjy];
ndoms=max(adj(:));

% Domain neighbors [E,W,N,S]
topo=zeros(ndoms,4);
topo(adjx(:,1),2)=adjx(:,2);
topo(adjx(:,2),1)=adjx(:,1);
topo(adjy(:,1),4)=adjy(:,2);
topo(adjy(:,2),3)=adjy(:,1);

% Domain interfaces [E,W,N,S]
net=zeros(ndoms,4);
net(adjx(:,1),2)=1:size(adjx,1);
net(adjx(:,2),1)=1:size(adjx,1);
net(adjy(:,1),4)=size(adjx,1)+(1:size(adjy,1));
net(adjy(:,2),3)=size(adjx,1)+(1:size(adjy,1));

% Horizontal interface neighbors [E,W]
RL=zeros(size(adjx));
[A,B]=ndgrid(adjx(:,1), adjx(:,2));
[I,J]=find(A==B);
RL(:,1)=1:size(adjx,1);
RL(I,2)=J;
RL(J,3)=I;

% Vertical interface neighbors [N,S]
TB=zeros(size(adjy));
[A,B]=ndgrid(adjy(:,1), adjy(:,2));
[I,J]=find(A==B);
TB(:,1)=(1:size(adjy,1))+size(adjx,1);
TB(I,2)=J+size(adjx,1);
TB(J,3)=I+size(adjx,1);

topo(topo==0)=ndoms+1;
net(net==0)=length(adj)+1;
RL(RL==0)=length(adj)+1;
TB(TB==0)=length(adj)+1;
end