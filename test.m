clear;

load test_cow.mat;



% d_max = max(dataPts,[],1);
% d_min = min(dataPts,[],1);
% L = d_max - d_min;
% bandwidth = max(L)/5;
% bandwidth = 0.75;

[clustCent,data2cluster,cluster2data,numClust] = MeanShiftCluster(dataPts,1,1.271111);

%  save('C:\\Users\\WYZ95\\desktop\\cluster.mat','data2cluster')