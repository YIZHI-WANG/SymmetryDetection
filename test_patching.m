load vertex.mat;
load face.mat;
load transformation.mat;
load vertex_in_cluster.mat;
load clusterSize.mat;

tic
[ patch1, patch2, patchSize] = Patching_entrance(5,0.05, vertex, face, transformation, vertex_in_cluster,clusterSize );
toc