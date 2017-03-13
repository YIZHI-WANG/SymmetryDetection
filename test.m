clear;

load test_cow_t.mat;
% testData = dataPts(1:2,:);
% save('C:\\Users\\WYZ95\\desktop\\testData.mat','testData')

% d_max = max(dataPts,[],1);
% d_min = min(dataPts,[],1);
% L = d_max - d_min;
% bandwidth = max(L)/5;
% bandwidth = 0.75;

[clustCent,data2cluster,cluster2data,numClust,clusterSize] = MeanShiftClusterT(dataPts,1,1.271111);

 save('C:\\Users\\WYZ95\\desktop\\cow_cluster0.mat','data2cluster')