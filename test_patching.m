load vertex.mat;
load face.mat;
load transformation.mat;
load vertex_in_cluster.mat;
load clusterSize.mat;


[ patch_count, patch1, patch2, patch_in_cluster_num, patch_2_clsuter,clsuter_2_patch,patchSize] = Patching_entrance(4,0.05, vertex, face, transformation, vertex_in_cluster,clusterSize );
