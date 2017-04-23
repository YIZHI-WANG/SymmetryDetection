function [ patch ] = Patching_core(threshold, patch, v1, v2, mode, scale, translation, vertex, face, v_ring, vf_ring)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



newMember = [];
%%find one-ring neighbors of v1
v1_vring = cell2mat(v_ring(v1));
[~, v1_vring_num] = size(v1_vring);
%%find unvisited neighbors
v1_vring_pos = vertex(:,v1_vring);
%%apply transformation
translated_v1_vring_pos = zeros(3,v1_vring_num);

if mode == 0
    translated_v1_vring_pos = v1_vring_pos.*scale + repmat(translation',[1 v1_vring_num]);
else
    translated_v1_vring_pos = (v1_vring_pos - repmat(translation',[1 v1_vring_num])) ./ scale;
end
%%find neighbor surface of v2
v2_vfring = cell2mat(vf_ring(v2));
[~, v2_vfring_num] = size(v2_vfring);
%%calculate dist(which dist)
dist = zeros(v2_vfring_num,v1_vring_num);
for v = 1:v1_vring_num
    v_pos = translated_v1_vring_pos(:,v);
    for f = 1:v2_vfring_num
        dist_ = compute_dist_of_point2tri( vertex, face, v_pos, v2_vfring(f) );
        dist(f,v) = dist_;
    end 
end

%%find min dist
dist = min(dist);
if mode == 0
    goodNeighbour = find(dist<=threshold*scale);
else
    goodNeighbour = find(dist<=threshold/scale);
end
% goodNeighbour = find(dist<=threshold);
%%add
newMember = setdiff(v1_vring(goodNeighbour),patch)';
patch = union(newMember, patch);
[row,col] = size(patch);
if col > row
    patch = patch';
end
[num,temp] = size(newMember);
if num > 0 && temp > 0
    for i=1:num
        newMember(i);
        patch = Patching_core(threshold, patch, newMember(i), v2, mode, scale, translation, vertex, face, v_ring, vf_ring);
    end
end
%%next ring

end

