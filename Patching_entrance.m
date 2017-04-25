function [ patch_count, patch1, patch2, patch_in_cluster_num, patch_2_clsuter,clsuter_2_patch,patchSize] = Patching_entrance( mode, threshold, vertex, face, transformation, vertex_in_cluster,clusterSize )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

save('D:\\vertex.mat','vertex')
save('D:\\face.mat','face')
save('D:\\transformation.mat','transformation')
save('D:\\vertex_in_cluster.mat','vertex_in_cluster')
save('D:\\clusterSize.mat','clusterSize')
    % threshold = norm(vertex(:,v1)+translation'-vertex(:,v2)) * 1.1;
%     poolSize = 8;N
%     if ischar(poolSize)  %where poolSize stores the user input
%     poolSize=str2num(poolSize); %convert string to number
%     end
%     parpool('local',poolSize);
    significant_mode = min(20,size(clusterSize,1));
%     significant_mode = size(clusterSize,1);
    boxSize = BoundingBoxSize(vertex);
    dist_threshold = boxSize * threshold;
    v_ring = compute_vertex_ring(face);
    vf_ring = compute_vertex_face_ring(face);
   
    
    [cluster_num,~] = size(vertex_in_cluster);
    [~,vertex_num] = size(vertex);
%     cluster_num = 1;
    
    
    
    P_patch1 = cell(significant_mode,10);
    P_patch2 = cell(significant_mode,10);
%     ones(significant_mode,10, vertex_num/2);
%     P_patch2 = ones(significant_mode,10, vertex_num/2);
%     P_patch1 = -1 .* P_patch1;
%     P_patch2 = -1 .* P_patch2;
    P_patchSize = cell(significant_mode,10);
    
    P_patch_num = zeros(1,significant_mode);
    
    [val,pos] = sort(clusterSize,'descend');
    height_threshold = val(significant_mode);
    
    total_patch_num = 0;
    
%     patching_cluster = sort(pos(1:20));
%     patching_cluster = patching_cluster';
%     parfor which_cluster = patching_cluster
tic
   parfor k = 1:significant_mode
        scale  = 1;
        which_cluster = pos(k);
        thisPatch1 = [];
        thisPatch2 = [];
        is_reflection = false;
        if clusterSize(which_cluster) >= height_threshold
            switch mode 
                case 0
                    translation = transformation(which_cluster,5:7);
                case 1
                    translation = transformation(which_cluster,4:6);
                case 3 
                    plane = transformation(which_cluster,1:4);
                case 4
                    rot = transformation(which_cluster,1:3);
                    plane = transformation(which_cluster,4:7);
                    is_reflection = reflection(rot);
                    if ~is_reflection
                         translation = plane(1:3) .* plane(4);
                    end
                case 5
                    scale = transformation(which_cluster,1);
                    rot = transformation(which_cluster,2:4);
                    plane = transformation(which_cluster,5:8);
                    is_reflection = reflection(rot);
                    if ~is_reflection
                         translation = plane(1:3) .* plane(4);
                    end
            end
            vertex_in_this_cluster = find(vertex_in_cluster(which_cluster,:)>0);
            vertex_in_this_cluster = vertex_in_cluster(which_cluster,vertex_in_this_cluster);
            [~,cluster_size] = size(vertex_in_this_cluster);
            pair_num = floor(cluster_size/2);
            
            patch1_in_cluster = -1 .* ones(pair_num, ceil(vertex_num/2));
            patch2_in_cluster = -1 .* ones(pair_num, ceil(vertex_num/2));
            patch_in_cluster = 0;
            
            for i = 1:pair_num
                confirmed = 1;
                v1 = vertex_in_this_cluster(2*i-1) + 1;
                v2 = vertex_in_this_cluster(2*i) + 1;
                skip = false;
%                 for p=1:patch_in_cluster
%                     v1_pos = find(patch1_in_cluster(p,:) == v1);
%                     v2_pos = find(patch2_in_cluster(p,:) == v2);
%                     if size(v1_pos,2) ==1 || size(v2_pos,2) == 1
%                         skip = true;
%                         break;
%                     end
%                 end
%                 
%                 if skip
%                     continue;
%                 end
                
                if mode  ==3 || mode == 4 || mode == 5
                    if is_reflection
                        dist_v1_2_plane = plane * [vertex(:,v1).*scale;1];
                        dist_v2_1_plane = plane * [vertex(:,v2);1];
                        if dist_v1_2_plane * dist_v2_1_plane >= 0
                            confirmed = 0;
                        end
                        translation = -2 * dist_v1_2_plane .*plane(1:3);   
                    end
                end
                v1_transformed = vertex(:,v1).*scale+ translation';
                dist = norm(v1_transformed-vertex(:,v2));
                
                if confirmed && dist <= dist_threshold*scale
                    %% get patches derived from v1 and v2
                    new_patch_1 = [];
                    new_patch_2 = [];
                    new_patch_1 = Patching_core(dist_threshold, new_patch_1, v1, v2, 0, scale, translation, vertex, face, v_ring, vf_ring);
                    new_patch_2 = Patching_core(dist_threshold, new_patch_2, v2, v1, 1, scale, translation, vertex, face, v_ring, vf_ring);
                    
                    
                    if size(new_patch_2,1) == 0 || size(new_patch_2,2) == 0 || size(new_patch_1,1) == 0 || size(new_patch_1,2) == 0
                        continue
                    end
                    %% merge pacthes if intersect with previous pacthes, which means spatial continunity
                    merge = false;
                    merge_with = 0;
                    for p=1:patch_in_cluster
                        common_vertex_1 = intersect(new_patch_1', patch1_in_cluster(p,:));
                        common_vertex_2 = intersect(new_patch_2', patch2_in_cluster(p,:));
                        if size(common_vertex_1,2)~= 0 && size(common_vertex_2,2)~= 0
                            merge_with = p;
                            merge = true;
                            break;
                        end
                    end
                    
                    if merge
                        union_patch1 = union(new_patch_1', patch1_in_cluster(merge_with,:));
                        union_patch2 = union(new_patch_2', patch2_in_cluster(merge_with,:));
                        if union_patch1(1) == -1
                            union_patch1 = union_patch1(2:size(union_patch1,2));
                        end
                        if union_patch2(1) == -1
                            union_patch2 = union_patch2(2:size(union_patch2,2));
                        end
                        patch1_in_cluster(p,1:size(union_patch1,2)) = union_patch1;
                        patch2_in_cluster(p,1:size(union_patch2,2)) = union_patch2;
                    else
                        patch_in_cluster = patch_in_cluster + 1;
                        patch1_in_cluster(patch_in_cluster,1:size(new_patch_1,1)) = new_patch_1';
                        patch2_in_cluster(patch_in_cluster,1:size(new_patch_2,1)) = new_patch_2';
                    end
                end
            end
            
            for p = 1:patch_in_cluster
                thisPatch1 = find(patch1_in_cluster(p,:)>0);
                thisPatch1 = patch1_in_cluster(p,thisPatch1);
                thisPatch2 = find(patch2_in_cluster(p,:)>0);
                thisPatch2 = patch2_in_cluster(p,thisPatch2);
                cut_patch1 = [];
                cut_patch2 = [];
                if is_reflection
                    vertex_in_patch1 = [vertex(:,thisPatch1); ones(1,size(thisPatch1,2))];
                    vertex_in_patch2 = [vertex(:,thisPatch2); ones(1,size(thisPatch2,2))];
                    dist_v1_2_plane = plane * vertex_in_patch1;
                    cut_patch1 = thisPatch1(find(dist_v1_2_plane<=0));
                    dist_v2_2_plane = plane * vertex_in_patch2;
                    cut_patch2 = thisPatch2(find(dist_v2_2_plane>=0));
                end
                thisPatch1 = setdiff(thisPatch1,cut_patch1);
                thisPatch2 = setdiff(thisPatch2,cut_patch2);
                temp1 = thisPatch1;
                temp2 = thisPatch2;
                thisPatch1 = setdiff(thisPatch1,temp2);
                thisPatch2 = setdiff(thisPatch2,temp1);
                thisPatch1Size = size(thisPatch1,2);
                thisPatch2Size = size(thisPatch2,2);
                if ceil(vertex_num/2) - thisPatch1Size > 0
                    complemetPatch1 = -1 .* ones(1,ceil(vertex_num/2) - thisPatch1Size);
                    newPatch1 = [thisPatch1 complemetPatch1];
                else
                    newPatch1 = thisPatch1;
                end

                if ceil(vertex_num/2) - thisPatch2Size > 0
                    complemetPatch2 = -1 .* ones(1,ceil(vertex_num/2) - thisPatch2Size);
                    newPatch2 = [thisPatch2 complemetPatch2];
                else
                    newPatch2 = thisPatch2;
                end
                if thisPatch1Size > 0 && thisPatch2Size > 0
                    P_patch1{k}{p} = newPatch1;
                    P_patch2{k}{p} = newPatch2;
                    P_patch_num(k) = P_patch_num(k) + 1;
                    P_patchSize{k}{p} = thisPatch1Size + thisPatch2Size;
                end
%                 
            end
        end
    end
toc
    patch_count = 0;
    for t=1:significant_mode
        patch_count = patch_count + P_patch_num(t);
    end  
    
    patch1 = zeros(patch_count,ceil(vertex_num/2));
    patch2 = zeros(patch_count,ceil(vertex_num/2));
    patch1 = -1 .* patch1;
    patch2 = -1 .* patch2;
    patch_2_clsuter = zeros(1,patch_count);
    clsuter_2_patch = zeros(cluster_num,10);
    patch_in_cluster_num = zeros(1,cluster_num);
    patchSize = zeros(1,patch_count);

    temp = 0;
    for m=1:significant_mode
        patch_in_cluster_num(pos(m)) = P_patch_num(m);
        for j=1:P_patch_num(m)
            temp = temp + 1;
            patch_2_clsuter(temp) = pos(m);
            clsuter_2_patch(pos(m),j) = temp;
            patch1(temp,:) = P_patch1{m}{j};
            patch2(temp,:) = P_patch2{m}{j};
            patchSize(temp) = P_patchSize{m}{j};
        end
    end  
    
    patch1 = patch1 - 1;
    patch2 = patch2 - 1;
    patch_2_clsuter = patch_2_clsuter -1;
    clsuter_2_patch = clsuter_2_patch -1;
%     parpool close
    end

function is_reflection = reflection(rot)
    if(rot(1) <3.1 || rot(1)>3.2)
        is_reflection = false;
        return
    end
    if(rot(2) <3.1 || rot(2)>3.2)
        is_reflection = false;
        return
    end
    if(rot(3) <3.1 || rot(3)>3.2)
        is_reflection = false;
        return
    end
    is_reflection = true;
end
