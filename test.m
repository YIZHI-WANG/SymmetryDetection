clear;

load baseballBat_RP7.mat;

[clustCent,data2cluster,cluster2dataMat,numClust,clusterSize,~] = MeanShiftCluster(dataPts,4,12,10);

% load dragon_T6.mat;
% 
% [clustCent,data2cluster,cluster2dataMat,numClust,clusterSize,~] = MeanShiftCluster(dataPts,1,1.27111,10);
% showResults(dataPts, data2cluster, 1);
%  save('C:\\Users\\WYZ95\\desktop\\cow_cluster0.mat','data2cluster')