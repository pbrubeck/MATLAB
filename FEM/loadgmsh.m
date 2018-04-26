function [tri, vert, bnd, edges] = loadgmsh(gmshfile)
enodes=zeros(51,1);
enodes([1:9,11,15])=[2;3;4;4;8;6;5;3;6;10;1];

fileID = fopen(gmshfile,'r');

fgetl(fileID); % $MeshFormat
fgetl(fileID); % X.Y filetype, data-size
fgetl(fileID); % $EndMeshFormat

fgetl(fileID); % $Nodes
nvert = fscanf(fileID, '%d', 1);
vert = fscanf(fileID, '%d %f %f %f\n', [4 nvert])';
vid = vert(:,1);
vert = vert(vid,2:3)*[1;1i];
fgetl(fileID); % $EndNodes

fgetl(fileID); % $Elements
ntri = fscanf(fileID, '%d', 1);
tri=zeros(ntri,6);
bnd=zeros(nvert,1);
edges=zeros(0,2);
for i=1:ntri
    eid = fscanf(fileID, '%d', 1);
    etype = fscanf(fileID, '%d', 1);
    ntags = fscanf(fileID, '%d', 1);
    tags  = fscanf(fileID, '%d', ntags)';
    k=enodes(etype);
    nodes = fscanf(fileID, '%d', k)';
    tri(eid,1:k) = nodes;
    fscanf(fileID, '\n', 1);
    
    if(etype==8)
        bnd(nodes)=1; 
        edges(end+1,:)=nodes([1,3]);
        edges(end+1,:)=nodes([2,3]);
    end
    
end

fgetl(fileID); % $EndElements
fclose(fileID);
bnd=find(bnd);
tri=tri(~any(tri==0,2),:);
end