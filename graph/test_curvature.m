% test for computation of curvature on meshes

% name = 'mushroom';
% name = 'fandisk';
% name = 'horse';
% name = 'nefertiti';
% name = 'elephant-50kv';
% name = 'bunny';
% name = 'skull';
% name = 'lion-head';
% name = 'vase-lion';
% name = 'david50kf';
% name = 'david_head';
% name = 'armadillo';
% name = 'aphro';
% name = 'gargoyle';
% name = 'tyra';
% name = 'screwdriver';
% name = 'cheburashka';
name = 'COW.OFF';
% name = 'remeshed.obj';
% options.name = name;
% 
% % path(path, '../off/');
% 
rep = ['results/curvature/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

[vertex,face] = read_mesh(name);
% [face_area,edge_length]=triangle(vertex,face);
% [ vertex, face ] = remesh( vertex, face ,mean(face_area)/1, 0 );
% write_obj("C:\\users\\WYZ95\\desktop\\remeshed.obj", vertex, face);

% load face.mat
% load vertex.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% local covariance analysis
options.covariance_smoothing = 15;
% [C,U,D] = compute_mesh_local_covariance(vertex,face,vertex,options);

% options for display
tau = 1.2;
options.normal_scaling = 1.5;


% options.normal = squeeze(U(:,2,:));
clf;
% options.face_vertex_color = perform_saturation( -D(2,:)' - D(3,:)',tau);
% plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-covariance.png'], 'png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% curvature

% size = BoundingBoxSize(vertex);
options.curvature_smoothing = 500;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);

options.normal = [];

figure(101)
clf;
options.face_vertex_color = perform_saturation(Cmax,tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cmax.png'], 'png');
figure(102)
clf;
options.face_vertex_color = perform_saturation(Cmin,tau);
plot_mesh(vertex,face, options);
shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cmin.png'], 'png');
% figure(3)
% clf;
% options.face_vertex_color = perform_saturation(Cmean,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% % saveas(gcf, [rep name '-cmean.png'], 'png');
% figure(4)
% clf;
% options.face_vertex_color = perform_saturation(Cgauss,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cgauss.png'], 'png');
% figure(5)
% clf;
% ratio = Cmin ./ Cmax;
% options.face_vertex_color = perform_saturation(ratio,tau);
% plot_mesh(vertex,face, options);
% shading interp; camlight; colormap jet(256);
% saveas(gcf, [rep name '-cabs.png'], 'png');
% dataPts = [Cmin Cmax]';
% [clustCent,data2cluster,cluster2dataMat,numClust,clusterSize] = MeanShiftCluster(dataPts,-1,0,10);


