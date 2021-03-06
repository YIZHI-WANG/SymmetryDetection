%testDistCluters

clear
profile on
% 
% nPtsPerClust = 2500;
% nClust  = 3;
% totalNumPts = nPtsPerClust*nClust;
% m(:,1) = [1 1]';
% m(:,2) = [-1 -1]';
% m(:,3) = [1 -1]';
% var = .6;
% bandwidth = 0.15;
% clustMed = [];
% %clustCent;
% 
% 
% x = var*randn(2,nPtsPerClust*nClust);
% %*** build the point set
% for i = 1:nClust
%     x(:,1+(i-1)*nPtsPerClust:(i)*nPtsPerClust)       = x(:,1+(i-1)*nPtsPerClust:(i)*nPtsPerClust) + repmat(m(:,i),1,nPtsPerClust);   
% end

load test_cow_t.mat;

% testData(1,:) = sum(abs(dataPts(1:3,:)));
% testData(2,:) = sum(abs(dataPts(4:6,:)));
% 
% testData(1,:) = testData(1,:)  ./ pi;
% testData(2,:) = testData(2,:) ./ 1.271111;


testData(1,:) = dataPts(1,:)./ pi ./ 2;
testData(2,:) = dataPts(4,:)./ 1.271111;


tic
[clustCent,point2cluster,~,~,~,clustMembsCell] = MeanShiftCluster(testData,2,0);
toc

numClust = length(clustMembsCell);


figure(10),clf,hold on
cVec = 'brcmykbrcmykbrcmykbrcmyk';%, cVec = [cVec cVec];
for k = 1:numClust-1
    myMembers = clustMembsCell{k};
    myClustCen = clustCent(:,k);
    plot(dataPts(1,myMembers),dataPts(4,myMembers),[cVec(mod(k,length(cVec)) + 1) '.'])
    plot(myClustCen(1) .*2.*pi  ,myClustCen(2) .*1.271111,'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(mod(k,length(cVec)) + 1), 'MarkerSize',10)
end
%%show noClusters
myMembers = clustMembsCell{numClust};
plot(dataPts(1,myMembers),dataPts(4,myMembers), ['g' 'o'])

title(['no shifting, numClust:' int2str(numClust)])