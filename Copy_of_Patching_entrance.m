function [ patch1, patch2, patchSize ] = Patching_entrance( mode, threshold, vertex, face, transformation, vertex_in_cluster,clusterSize )
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
    significant_mode = min(50,size(clusterSize,1));
%     significant_mode = size(clusterSize,1);
    boxSize = BoundingBoxSize(vertex);
    dist_threshold = boxSize * threshold;
    v_ring = compute_vertex_ring(face);
    vf_ring = compute_vertex_face_ring(face);
   
    
    [cluster_num,~] = size(vertex_in_cluster);
    [~,vertex_num] = size(vertex);
%     cluster_num = 1;
    
    
    
    P_patch1 = ones(significant_mode, vertex_num);
    P_patch2 = ones(significant_mode, vertex_num);
    P_patch1 = -1 .* P_patch1;
    P_patch2 = -1 .* P_patch2;
    P_patchSize = zeros(1,significant_mode);
    
    [val,pos] = sort(clusterSize,'descend');
    height_threshold = val(significant_mode);
    
    
    
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
        %             scale = transformation(which_cluster,1);
        %             rotation = transformation(which_cluster,2:4);
                    translation = transformation(which_cluster,5:7);
                case 1
        %             rotation = transformation(which_cluster,1:3);
                    translation = transformation(which_cluster,4:6);
        %             scale = 1;
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
            pair_num = cluster_size/2;
            for i = 1:pair_num
                confirmed = 1;
                v1 = vertex_in_this_cluster(2*i-1) + 1;
                v2 = vertex_in_this_cluster(2*i) + 1;
                if mode  ==3 || mode == 4 || mode == 5
                    if is_reflection
                        dist_v1_2_plane = plane * [vertex(:,v1).*scale;1];
                        dist_v2_1_plane = plane * [vertex(:,v2);1];
                        if dist_v1_2_plane * dist_v2_1_plane >= 0
                            confirmed = 0;
                        end
                        translation = -2 * dist_v1_2_plane .*plane(1:3);   
%                     if dist_v1_2_plane > 0
%                         translation = -2 * dist_v1_2_plane .*plane(1:3);
%                     else
%                         translation = 2 * dist_v1_2_plane .*plane(1:3);
%                     end
                    end

                end
                v1_transformed = vertex(:,v1).*scale+ translation';
                dist = norm(v1_transformed-vertex(:,v2));
%                 dist_threshold = norm(translation) * threshold;
                if confirmed && dist <= dist_threshold*scale
%                 if confirmed && dist <= dist_threshold
                    if find(thisPatch1 == v1)
                        %%v1 has been in this patch
    %                     test = 1
                    else
%                         thisPatch1 = [];
                        thisPatch1 = union(thisPatch1, v1);
                        thisPatch1 = Patching_core(dist_threshold, thisPatch1, v1, v2, 0, scale, translation, vertex, face, v_ring, vf_ring);
                    end
                    if find(thisPatch2 == v2)
                        %%v2 has been in this patch
    %                     test = 2
                    else
                        thisPatch2 = union(thisPatch2, v2);
                        thisPatch2 = Patching_core(dist_threshold, thisPatch2, v2, v1, 1, scale, translation, vertex, face, v_ring, vf_ring);
                    end 
                end
            end
        end
        
        [thisPatch1Size,~] = size(thisPatch1);
        [thisPatch2Size,~] = size(thisPatch2);
        
        cut_patch1 = [];
        cut_patch2 = [];
        if is_reflection
            for m = 1:thisPatch1Size
                this_vertex = thisPatch1(m);
                dist_2_plane = plane * [vertex(:,this_vertex).*scale;1];
                if dist_2_plane <=0
                    cut_patch1 = union(cut_patch1,this_vertex);
                end
            end
            for m = 1:thisPatch2Size
                this_vertex = thisPatch2(m);
                dist_2_plane = plane * [vertex(:,this_vertex);1];
                if dist_2_plane >=0
                    cut_patch2 = union(cut_patch2,this_vertex);
                end
            end
            thisPatch1 = setdiff(thisPatch1,cut_patch1);
            thisPatch2 = setdiff(thisPatch2,cut_patch2);
        end
        
        temp1 = thisPatch1;
        temp2 = thisPatch2;
        thisPatch1 = setdiff(thisPatch1,temp2);
        thisPatch2 = setdiff(thisPatch2,temp1);
        [thisPatch1Size,~] = size(thisPatch1);
        [thisPatch2Size,~] = size(thisPatch2);
%         which_cluster
%         "start"
        if vertex_num - thisPatch1Size > 0
            complemetPatch1 = -1 .* ones(1,vertex_num - thisPatch1Size);
            newPatch1 = [thisPatch1' complemetPatch1];
        else
            newPatch1 = thisPatch1';
        end
        
        if vertex_num - thisPatch2Size > 0
            complemetPatch2 = -1 .* ones(1,vertex_num - thisPatch2Size);
            newPatch2 = [thisPatch2' complemetPatch2];
        else
            newPatch2 = thisPatch2';
        end
%         if thisPatch1Size > 0 && thisPatch2Size > 0
        if thisPatch1Size > 0 
            P_patch1(k,:) = newPatch1;
        end
        if thisPatch2Size > 0 
            P_patch2(k,:) = newPatch2;
        end
        P_patchSize(1,k) = thisPatch1Size + thisPatch2Size;
%         end
%         which_cluster
%         "end"
    end
toc
    patch1 = ones(cluster_num,vertex_num);
    patch2 = ones(cluster_num,vertex_num);
    patch1 = -1 .* patch1;
    patch2 = -1 .* patch2;
    patchSize = zeros(1,cluster_num);
    
    for k=1:significant_mode
        patch1(pos(k),:) = P_patch1(k,:);
        patch2(pos(k),:) = P_patch2(k,:);
        patchSize(1,pos(k)) = P_patchSize(1,k);
    end  
    
    patch1 = patch1 - 1;
    patch2 = patch2 - 1;
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
