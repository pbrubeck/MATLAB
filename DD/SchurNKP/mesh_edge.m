classdef mesh_edge
    %MESH_EDGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        face1
        face2
    end
    
    methods
        function obj = mesh_edge(f1,f2)
            if nargin>0
                obj.face1=f1;
            else
                obj.face1=[-1,-1];
            end
            if nargin>1
                obj.face2=f2;
            else
                obj.face2=[-1,-1];
            end
        end
        
        function obj = add_face(obj,face)
            if all(obj.face1 == [-1,1])
                obj.face1=face;
            else
                obj.face2=face;
            end
        end
        
        function b=initialized(obj)
             b = all(obj.face1 ~= [-1,-1]) && all(obj.face2 ~= [-1,-1]);
        end
        
        function h=hash(obj)
            h = obj.face1 + obj.face2;
        end
        
    end
end

