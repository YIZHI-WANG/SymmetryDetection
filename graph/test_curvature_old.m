
load face.mat
load vertex.mat

%% curvature

options.curvature_smoothing = BoundingBoxSize(vertex);
options.null = 0;
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);
