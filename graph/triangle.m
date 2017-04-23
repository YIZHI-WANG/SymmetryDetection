function [face_area,edge_length]=triangle(vertex,face)
    face_num = size(face,2);
    face_area = zeros(face_num,1);
    edge_length = zeros(face_num,3);
    for i = 1:face_num
        v1 = vertex(:,face(1,i));
        v2 = vertex(:,face(2,i));
        v3 = vertex(:,face(3,i));
        % Length of edges
        L=[sqrt(sum((v1-v2).^2)) sqrt(sum((v2-v3).^2)) sqrt(sum((v3-v1).^2))];

        % Area calculation with Heron's formula
        s = ((L(1)+L(2)+L(3))/2); 
        face_area(i) = sqrt(s*(s-L(1))*(s-L(2))*(s-L(3))); 
        edge_length(i,:) = L;
    end
end 