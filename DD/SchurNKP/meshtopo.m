function [net, adj, corners, edges, bnd] = meshtopo(quads)
nquad=size(quads,1);
[adj,bnd]=faces_to_edges(quads);
[net,adj]=adjacency(adj, nquad);

[verts_faces, verts_vertices]=verts_from_quad(quads);
loops=get_loops(verts_faces, verts_vertices);
corners=loops_to_corners(loops, quads);
edges=corners_to_edges(adj, corners);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        Edge Functionality
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edge_key]=face_to_key(face)
    edge_key = sprintf('%d,%d', [min(face), max(face)]);
end

function edge=add_face_to_edge(edge_map, key, face)
    if edge_map.isKey(key)
        edge = edge_map(key);
        edge(3:4) = face;
    else
        edge = [face, 0, 0];
    end
end

function [interfaces, boundaries]=faces_to_edges(quad)
    edge_map = containers.Map;
    side = [1,3; 2,4; 1,2; 3,4];
    for q = 1:size(quad,1)
        for s = 1:size(side,1)
            key = face_to_key(quad(q, side(s,:)));
            edge_map(key) = add_face_to_edge(edge_map, key, [s,q]);
        end
    end
    edges_list = cell2mat(reshape(edge_map.values(),[],1));
    p = any(edges_list==0, 2);
    interfaces = edges_list(~p,:);
    boundaries = edges_list(p,1:2);
end

function [net,adj]=adjacency(adj, nquad)
% This sorts the interface numbering such that we get minimum fill in the
% Schur complement.
[x11,y11]=ndgrid(adj(:,2), adj(:,2));
[x12,y12]=ndgrid(adj(:,2), adj(:,4));
[x21,y21]=ndgrid(adj(:,4), adj(:,2));
[x22,y22]=ndgrid(adj(:,4), adj(:,4));
mask=(x11==y11)+(x12==y12)+(x21==y21)+(x22==y22);
p=symamd(mask);
adj=adj(p,:);

s1=adj(:,1);
q1=adj(:,2);
s2=adj(:,3);
q2=adj(:,4);
net=zeros(nquad,4);
for s=1:4
    net(q1(s1==s),s)=find(s1==s);
    net(q2(s2==s),s)=find(s2==s);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        Vertex Functionality
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b]=is_vertex_loop(vert_vertices, verts_faces)
    b= numel(vert_vertices)==numel(verts_faces);
end

function [loop_verts_nums]=get_loops(verts_faces, verts_vertices)
    loop_verts_nums=zeros(numel(verts_faces),1);
    num_loops = 0;
    for vert_id=1:numel(verts_faces)
        if is_vertex_loop(verts_faces{vert_id}, verts_vertices{vert_id})
            num_loops = num_loops+1;
            loop_verts_nums(num_loops)=vert_id;
        end
    end
    loop_verts_nums=loop_verts_nums(1:num_loops);
end

function [corners] = loops_to_corners(loops, quads)
    corners=zeros(size(quads));
    for k = 1:numel(loops)
        corners(quads==loops(k))=k;
    end
end

function [verts_faces] = add_face_to_vertex_faces(verts_faces, face)
    b = false;
    for i=1:numel(verts_faces)
        b = b || all(verts_faces{i} == face);
    end
    if (~b)
        verts_faces{end+1,1} = face;
    end
end

function [verts_vertices] = add_vertex_to_vertex_vertices(verts_vertices, vertex)
    b = false;
    for i=1:numel(verts_vertices)
        b = b || all(verts_vertices{i} == vertex);
    end
    if (~b)
        verts_vertices{end+1,1} = vertex;
    end
end

function [verts_faces, verts_vertices] = verts_from_quad(quads)
    nquads = size(quads,1);
    num_verts=0;
    verts_faces = cell(4*nquads,1);
    [verts_faces{:}] = deal(cell(0,1));
    
    verts_vertices = cell(4*nquads,1);
    [verts_vertices{:}] = deal(cell(0,1));
    
    for face_idx=1:nquads
        face_verts = quads(face_idx,:);
        for corner=1:numel(face_verts)
            vert = face_verts(corner);
            
            % Add attached face to the vertex
            verts_faces{vert} = add_face_to_vertex_faces(verts_faces{vert}, [corner, face_idx]);
            
            % Add the neighboring vertices to the vertix
            verts_vertices{vert}= add_vertex_to_vertex_vertices(verts_vertices{vert},  face_verts(bitxor(corner-1,1)+1));
            verts_vertices{vert}= add_vertex_to_vertex_vertices(verts_vertices{vert},  face_verts(bitxor(corner-1,2)+1));
            
            num_verts = max(num_verts, vert);
        end
    end
    
    verts_faces=verts_faces(1:num_verts);
    verts_vertices=verts_vertices(1:num_verts);
end

function edges=corners_to_edges(adj, corners)
    edges=zeros(size(adj,1),2);
    east =adj(:,1)==1;
    west =adj(:,1)==2;
    north=adj(:,1)==3;
    south=adj(:,1)==4;
    
    edges(east, :)=corners(adj(east ,2), [1,3]);
    edges(west, :)=corners(adj(west ,2), [2,4]);
    edges(north,:)=corners(adj(north,2), [1,2]);
    edges(south,:)=corners(adj(south,2), [3,4]);
end