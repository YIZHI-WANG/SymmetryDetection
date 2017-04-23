function [ new_vertexs, new_faces, new_vertex_num, new_face_num] = remeshing_core(vertex, face ,threshold)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% A = 
next = [2 3 1];

[face_area,edge_length] = triangle(vertex,face);
edges = compute_edges(face);
e2f = compute_edge_face_ring(face);
cut_faces = find(face_area >= threshold);
cut_num = size(cut_faces,1);
edge_num = size(edges,2); %%TODO edge num
vertex_num = size(vertex,2);
face_num = size(face,2);
new_vertex_num = 0;
new_face_num = 0;
cut_dictionary = sparse(vertex_num,vertex_num);
face_cut_sign = zeros(face_num,1);

for i = 1:cut_num
    cut_face = cut_faces(i);
%     [~,longest_edge_num] = max(edge_length(i,:));
    for edge = 1:3
        v1 = face(edge,cut_face);
        v2 = face(next(edge),cut_face);
        if v2 < v1
            temp = v1;
            v1 = v2;
            v2 = temp;
        end
        if cut_dictionary(v1,v2) == 0
            new_vertex_num = new_vertex_num+1;
            cut_dictionary(v1,v2) = new_vertex_num+vertex_num;
    %         cut_dictionary(v2,v1) = new_vertex_num;
        end
    end
end

[cut_edge_v1,cut_edge_v2] = find(cut_dictionary ~= 0);

for i = 1:new_vertex_num
    v1 = cut_edge_v1(i);
    v2 = cut_edge_v2(i);
%     vertex_num = vertex_num + 1;
    which_vertex = max(cut_dictionary(v1,v2),cut_dictionary(v2,v1));
    new_vertex = (vertex(:,v1)+vertex(:,v2))./2;
    vertex(:,which_vertex) = new_vertex;
    face1_num = e2f(v1,v2);
    face2_num = e2f(v2,v1);
    if face1_num > 0
%         face1 = face(:,face1_num);
        face_cut_sign(face1_num) = 1;
%         v1_pos = find(face1 == v1);
%         v2_pos = find(face1 == v2);
%         face(v1_pos,face1_num) = vertex_num;
%         face_num = face_num+1;
%         face(:,face_num) = face1;
%         face(v2_pos,face_num) = vertex_num;
%         new_face_num = new_face_num + 1;
    end
    if face2_num > 0
%         face2 = face(:,face2_num);
        face_cut_sign(face2_num) = 1;
%         v1_pos = find(face2 == v1);
%         v2_pos = find(face2 == v2);
%         face(v1_pos,face2_num) = vertex_num;
%         face_num = face_num+1;
%         face(:,face_num) = face2;
%         face(v2_pos,face_num) = vertex_num;
%         new_face_num = new_face_num + 1;
    end
end

new_face = zeros(3,2*new_vertex_num); 
need_cut_faces = find(face_cut_sign == 1);
need_cut_faces_num = size(need_cut_faces,1);
for j = 1:need_cut_faces_num
    cut_face = face(:,need_cut_faces(j));
    face_v1 = cut_face(1);
    face_v2 = cut_face(2);
    face_v3 = cut_face(3);
    cut_point = [0 0 0];
    if cut_dictionary(face_v1,face_v2) ~= 0 || cut_dictionary(face_v2,face_v1) ~= 0
        cut_point(1) = 1;
    end
    if cut_dictionary(face_v2,face_v3) ~= 0 || cut_dictionary(face_v3,face_v2) ~= 0
        cut_point(2) = 1;
    end
    if cut_dictionary(face_v3,face_v1) ~= 0 || cut_dictionary(face_v1,face_v3) ~= 0
        cut_point(3) = 1;
    end
    [cut_point,cut_point_pos] = find(cut_point==1);
    switch size(cut_point,2)
        case 1
            cut_v1 = cut_face(cut_point_pos);
            cut_v2 = cut_face(next(cut_point_pos));
            mid_point = max(cut_dictionary(cut_v1,cut_v2),cut_dictionary(cut_v2,cut_v1));
            new_face_num = new_face_num +1;
            new_face(:,new_face_num) = cut_face;
            new_face(cut_point_pos,new_face_num) = mid_point;
            new_face_num = new_face_num +1;
            new_face(:,new_face_num) = cut_face;
            new_face(next(cut_point_pos),new_face_num) = mid_point;
        case 2
            no_cut_pos = setdiff([1 2 3], cut_point_pos);
            switch no_cut_pos
                case 1
                    mid_point_1 = max(cut_dictionary(face_v2,face_v3),cut_dictionary(face_v3,face_v2));
                    mid_point_2 = max(cut_dictionary(face_v3,face_v1),cut_dictionary(face_v1,face_v3));
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [mid_point_1 face_v3 mid_point_2]';
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [mid_point_1 mid_point_2 face_v2]';
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [face_v1 face_v2 mid_point_2]';
                case 2
                    mid_point_1 = max(cut_dictionary(face_v3,face_v1),cut_dictionary(face_v1,face_v3));
                    mid_point_2 = max(cut_dictionary(face_v1,face_v2),cut_dictionary(face_v2,face_v1));
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [mid_point_1 face_v1 mid_point_2]';
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [mid_point_1 mid_point_2 face_v3]';
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [face_v2 face_v3 mid_point_2]';
                case 3
                    mid_point_1 = max(cut_dictionary(face_v1,face_v2),cut_dictionary(face_v2,face_v1));
                    mid_point_2 = max(cut_dictionary(face_v2,face_v3),cut_dictionary(face_v3,face_v2));
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [mid_point_1 face_v2 mid_point_2]';
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [mid_point_1 mid_point_2 face_v1]';
                    new_face_num = new_face_num +1;
                    new_face(:,new_face_num) = [face_v3 face_v1 mid_point_2]';
            end
            
        case 3
            mid_point_1 = max(cut_dictionary(face_v1,face_v2),cut_dictionary(face_v2,face_v1));
            mid_point_2 = max(cut_dictionary(face_v2,face_v3),cut_dictionary(face_v3,face_v2));
            mid_point_3 = max(cut_dictionary(face_v3,face_v1),cut_dictionary(face_v1,face_v3));
            new_face_num = new_face_num +1;
            new_face(:,new_face_num) = [face_v1 mid_point_1 mid_point_3]';
            new_face_num = new_face_num +1;
            new_face(:,new_face_num) = [mid_point_1 face_v2 mid_point_2]';
            new_face_num = new_face_num +1;
            new_face(:,new_face_num) = [mid_point_2 face_v3 mid_point_3]';
            new_face_num = new_face_num +1;
            new_face(:,new_face_num) = [mid_point_2 mid_point_3 mid_point_1]';
    end
    
end

    new_vertexs = vertex;
    reserve_faces = find(face_cut_sign == 0);
    new_faces = [face(:,reserve_faces) new_face(:,1:new_face_num)];
% new_faces = [face(:,reserve_faces)];

end







