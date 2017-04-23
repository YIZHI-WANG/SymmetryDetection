function [ new_vertexs, new_faces ] = remesh( vertex, face ,threshold, mode )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if mode == 1
    face = face';
    vertex = vertex';
end
[new_vertexs,new_faces,new_vertex_num,new_face_num] = remeshing_core(vertex, face ,threshold);

while new_face_num > 0
    [ new_vertexs, new_faces, new_vertex_num, new_face_num] = remeshing_core(new_vertexs, new_faces,threshold);
end

write_obj("C:\\users\\yizhiw\\desktop\\remeshed.obj", new_vertexs', new_faces');

if mode == 1
    new_faces = new_faces';
    new_vertexs = new_vertexs';
end
% new_faces = face;
% new_vertexs = vertex;


end

