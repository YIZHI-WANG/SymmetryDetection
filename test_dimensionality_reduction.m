clc
clear
close all
 
load test_cow_t.mat;
load cow_cluster0.mat;

boxSize = 1.08571;
metric = dataPts';
% metric(:,1:3)= metric(:,1:3) .* (boxSize / 2 / pi);

% PCA降维
[mappedX, mapping] = compute_mapping(metric, 'PCA', 2);
x_max = max(mappedX(:,1));
x_min = min(mappedX(:,1));
y_max = max(mappedX(:,2));
y_min = min(mappedX(:,2));
clusterSize = max(data2cluster);
color = jet(clusterSize);
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';
figure
for i = 1:clusterSize
    thisCluster = find(data2cluster == i);
    
    plot(mappedX(thisCluster,1), mappedX(thisCluster,2),'o','markerfacecolor',cVec(mod(i,length(cVec)) + 1))
    axis([x_min x_max y_min  y_max]);
    title(['Cluster '  num2str(i)])
%     hold on
    pause
%     clf
end

% figure
% scatter(mappedX(:,1), mappedX(:,2),3,data2cluster,'filled')
% title('Result of PCA')


% % 产生测试数据
% [X, labels] = generate_data('helix', 2000);
% figure
% scatter3(X(:,1), X(:,2), X(:,3), 5, labels)
% title('Original dataset')
% drawnow
%  
% % 估计本质维数
% no_dims = round(intrinsic_dim(X, 'MLE'));
% disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
 
% % PCA降维
% [mappedX, mapping] = compute_mapping(X, 'PCA', no_dims);
% figure
% scatter(mappedX(:,1), mappedX(:,2), 5, labels)
% title('Result of PCA')
%  
% % Laplacian降维
% [mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
% figure
% scatter(mappedX(:,1), mappedX(:,2), 5, labels(mapping.conn_comp))
% title('Result of Laplacian Eigenmaps')
% drawnow
%  
% % % Isomap降维
% % [mappedX, mapping] = compute_mapping(X, 'Isomap', no_dims);
% % figure
% % scatter(mappedX(:,1), mappedX(:,2), 5, labels(mapping.conn_comp))
% % title('Result of Isomap')
% % drawnow