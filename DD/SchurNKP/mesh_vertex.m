classdef mesh_vertex
    %MESH_VERTEX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        faces
        vertices
        num
    end
    
    methods
        function obj = mesh_vertex(num, faces, vertices)
            obj.num=num;
            obj.faces=faces;
            obj.vertices=vertices;
        end
        
        function obj = add_face(obj,face)
           obj.faces(end+1)=face;
        end
        
        function obj = add_vertex(obj,vertex)
           obj.vertices(end+1)=vertex;
        end
        
        function b = is_loop(obj)
           b=size(obj.vertices,1)==size(obj.faces,1);
        end
    end
end

