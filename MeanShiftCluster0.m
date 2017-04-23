function [clustCent,data2cluster,cluster2dataMat,numClust,clusterSize,cluster2dataCell] = MeanShiftCluster(dataPts,mode,boxSize,grid)
tic

zeroCluster = false;
strictVote = false;

switch mode
    case 0
        save('D:\\dragon_T7_0.mat','dataPts')
%         maxScale = max(dataPts(1,:));
%         maxRot = max(max(dataPts(2:4,:)) -  min(dataPts(2:4,:)));
%         maxTrans = max(max(dataPts(5:7,:)) -  min(dataPts(5:7,:)));
%         maxRotx =  max(max(dataPts(2,:)) -  min(dataPts(2,:)));
%         maxRoty =  max(max(dataPts(4,:)) -  min(dataPts(4,:)));
%         maxRotz =  max(max(dataPts(4,:)) -  min(dataPts(4,:)));
%         maxTransx =  max(max(dataPts(5,:)) -  min(dataPts(5,:)));
%         maxTransy =  max(max(dataPts(6:)) -  min(dataPts(6,:)));
% %         maxTransz =  max(max(dataPts(7,:)) -  min(dataPts(7,:)));
%         dataPts(1,:) =dataPts(1,:) ./ maxScale;
% %         dataPts(2:4,:) = dataPts(1:3,:)  ./ pi ./ 2;
%         dataPts(2:4,:) = dataPts(2:4,:)  ./maxRot;
%         dataPts(5:7,:) = dataPts(5:7,:)  ./ maxTrans ;
        dataPts(1,:) =dataPts(1,:) ./ 5;
        dataPts(2:4,:) = dataPts(2:4,:)  ./ pi ./ 2;
        dataPts(5:7,:) = dataPts(5:7,:)  ./ boxSize ;
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L(1,2:7))/grid;
    case 1
        save('D:\\dragon_T6.mat','dataPts')
        dataPts(1:3,:) = dataPts(1:3,:)  ./ pi ./ 2;
        dataPts(4:6,:) = dataPts(4:6,:)  ./ boxSize ;
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L)/grid;
    case 2
        save('D:\\baseballBat_R3.mat','dataPts')
        dataPts(1:3,:) = dataPts(1:3,:)  ./ pi ./ 2;
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L)/grid;
    case 3
        save('D:\\baseballBat_P4.mat','dataPts')
        dataPts(4,:) = dataPts(4,:)  ./ boxSize;
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L)/grid;
    case 4
        save('D:\\baseballBat_RP7.mat','dataPts')
        dataPts(1:3,:) = dataPts(1:3,:)  ./ pi ./ 2;
        dataPts(7,:) = dataPts(7,:)  ./ boxSize;
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L)/grid;
    case 5
        save('D:\\baseballBat_SV3.mat','dataPts')
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L)/grid;
    case 6
        save('D:\\baseballBat_RV6.mat','dataPts')
        dataPts(1:3,:) = dataPts(1:3,:)  ./ pi ./ 2;
        d_max = max(dataPts',[],1);
        d_min = min(dataPts',[],1);
        L = d_max - d_min;
        bandWidth = max(L)/grid;
    case 7
        bandWidth = 0.1;
end

plotFlag = false;

%**** Initialize stuff ***
[numDim,numPts] = size(dataPts);
numClust        = 0;
bandSq          = bandWidth.^2;
initPtInds      = 1:numPts;
maxPos          = max(dataPts,[],2);                          %biggest size in each dimension
minPos          = min(dataPts,[],2);                          %smallest size in each dimension
boundBox        = maxPos-minPos;                        %bounding box size
sizeSpace       = norm(boundBox);                       %indicator of size of data space
stopThresh      = 1e-3*bandWidth;                       %when mean has converged
clustCent       = [];                                   %center of clust
beenVisitedFlag = zeros(1,numPts,'uint8');              %track if a points been seen already
numInitPts      = numPts;                               %number of points to posibaly use as initilization points
clusterVotes    = zeros(1,numPts,'uint16');             %used to resolve conflicts on cluster membership

% set these weights so that a rotation by 180 degrees corresponds to a
% displacement of half the bounding box diagonal
Beta1 = 10;
Beta2 = (boxSize / 2 / pi) ^ 2;
Beta3 = 1;

if zeroCluster
    myZeroMean          = zeros(numDim,1);                     % intilize mean to this points location                        
    metric = (repmat(myZeroMean,1,numPts) - dataPts).^2;
    sqDistToAll = sum(metric);    %dist squared from mean to all points still active
    zeroMembers      = find(sqDistToAll < bandSq);               %points within bandWidth
    beenVisitedFlag(zeroMembers) = 1;                         %mark that these points have been visited
    zeroClusterVotes    = zeros(1,numPts,'uint16');         % used to resolve conflicts on cluster membership
    zeroClusterVotes(zeroMembers) = zeroClusterVotes(zeroMembers)+1;  %add a vote for all the in points belonging to this cluster
    notZeroMembers = setdiff(1:numPts, zeroMembers);
    initPtInds = notZeroMembers;
    numInitPts = size(notZeroMembers,2);
end

while numInitPts

    tempInd         = ceil( (numInitPts-1e-6)*rand);        % pick a random seed point
    stInd           = initPtInds(tempInd);                  % use this point as start of mean
    myMean          = dataPts(:,stInd);                     % intilize mean to this points location
    myMembers       = [];                                   % points that will get added to this cluster                          
    thisClusterVotes    = zeros(1,numPts,'uint16');         % used to resolve conflicts on cluster membership

    while 1     %loop untill convergence
        metric = (repmat(myMean,1,numPts) - dataPts).^2;
        sqDistToAll = sum(metric);    %dist squared from mean to all points still active
        inInds      = find(sqDistToAll < bandSq);               %points within bandWidth

        if zeroCluster
        inInds = setdiff(inInds, zeroMembers);
        end
        
        thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;  %add a vote for all the in points belonging to this cluster
        
        myOldMean   = myMean;                                   %save the old mean
        
        kernel_value = EpanechnikovKernel(numDim, size(inInds,2), sqDistToAll(inInds)./bandSq);
        numerator  = sum(dataPts(:,inInds) * kernel_value' , 2);
        denominator = sum(kernel_value , 2)+0.0000000001;
        if isnan(myMean(1,1))
             disp('==NaN!')
        end
        if denominator == 0
            disp('==0!')
        end
        
        if max(isnan(myMean))
            disp('here')
        end
        
        myMean = numerator ./ denominator;

%         myMean      = mean(dataPts(:,inInds),2);                %compute the new mean
       
        myMembers   = [myMembers inInds];                       %add any point within bandWidth to the cluster
        beenVisitedFlag(myMembers) = 1;                         %mark that these points have been visited
        
        %*** plot stuff ****
        if plotFlag
            figure(12345),clf,hold on
            if numDim == 2
                plot(dataPts(1,:),dataPts(2,:),'.')
                plot(dataPts(1,myMembers),dataPts(2,myMembers),'ys')
                plot(myMean(1),myMean(2),'go')
                plot(myOldMean(1),myOldMean(2),'rd')
                movement = norm(myMean - myOldMean)
                ratio = movement / bandWidth
                plot([myOldMean(1) myMean(1)],[myOldMean(2) myMean(2)],'b-')
                pause
            end
        end

        %**** if mean doesn't move much stop this cluster ***
%         test = norm(myMean-myOldMean)
%         if isnan(test)
%             disp('fail')
%         end
        dist = norm(myMean-myOldMean);
        if dist < stopThresh
            
            %check for merge posibilities
            mergeWith = 0;
            for cN = 1:numClust
                distToOther = norm(myMean-clustCent(:,cN));     %distance from posible new clust max to old clust max
                if distToOther < bandWidth/2                    %if its within bandwidth/2 merge new and old
                    mergeWith = cN;
%                     mergeWith = 0;
                    break;
                end
            end
            
            
            if mergeWith > 0    % something to merge
%                 clustCent(:,mergeWith)       = 0.5*(myMean+clustCent(:,mergeWith));             %record the max as the mean of the two merged (I know biased twoards new ones)
                
                myMember = find(thisClusterVotes>0);
                [~,myMemberSize] = size(myMember);
                hisMember = find( clusterVotes(mergeWith,:)>0);
                [~,hisMemberSize] = size(hisMember);
                clustCent(:,mergeWith) = (myMemberSize .* myMean + hisMemberSize .*clustCent(:,mergeWith) ) ./ (myMemberSize + hisMemberSize);
                %clustMembsCell{mergeWith}    = unique([clustMembsCell{mergeWith} myMembers]);   %record which points inside 
                clusterVotes(mergeWith,:)    = clusterVotes(mergeWith,:) + thisClusterVotes;    %add these votes to the merged cluster
            else    %its a new cluster
                numClust                    = numClust+1;                   %increment clusters
                clustCent(:,numClust)       = myMean;                       %record the mean  
                %clustMembsCell{numClust}    = myMembers;                    %store my members
                clusterVotes(numClust,:)    = thisClusterVotes;
            end

            break;
        end

    end
    
    
    initPtInds      = find(beenVisitedFlag == 0);           %we can initialize with any of the points not yet visited
    numInitPts      = length(initPtInds);                   %number of active points in set

end

if strictVote
    for i = 1:numClust
        voters = find(clusterVotes(i,:) > 0);
        meanVote = mean(clusterVotes(i,voters));
        badVoters = find(clusterVotes(i,:) < 0.5 * meanVote);
        clusterVotes(i,badVoters) = 0;
    end
end


if zeroCluster
    numClust                    = numClust+1;                   %increment clusters
    clustCent(:,numClust)       = myZeroMean;                       %record the mean  
    %clustMembsCell{numClust}    = myMembers;                    %store my members
    clusterVotes(numClust,:)    = inf *zeroClusterVotes ;
end


                
clusterSize = zeros(numClust,1);
[dataVotes,data2cluster] = max(clusterVotes,[],1);                %a point belongs to the cluster with the most votes


%%find no cluster data
if strictVote
    badData = find(dataVotes == 0);
    data2cluster(badData) = -1;
end

cluster2dataMat = -1 .* ones(numClust,numPts);
cluster2dataCell = cell(numClust,1);


if mode ~= 7
%*** If they want the cluster2data cell find it for them
    if nargout > 2
    %     cluster2dataCell = cell(numClust,1);
        for cN = 1:numClust
            myMembers = find(data2cluster == cN);
            [~,cluster_size] = size(myMembers);
            clusterSize(cN) =  cluster_size;
             cluster2dataMat(cN,1:cluster_size) = myMembers;
    %         cluster2dataCell{cN} = myMembers;
        end
    end
    cluster2dataMat = cluster2dataMat - 1;
    data2cluster = data2cluster - 1;
    data2cluster = data2cluster';
    
    switch mode
        case 0 
            clustCent(1,:) = clustCent(1,:) .* 5;
            clustCent(2:4,:) = clustCent(2:4,:) .* pi .* 2 ;
            clustCent(5:7,:) = clustCent(5:7,:) .* boxSize ;
%             clustCent(1,:) = clustCent(1,:) .* maxScale;
%             clustCent(2:4,:) = clustCent(2:4,:) .* maxRot ;
%             clustCent(5:7,:) = clustCent(5:7,:) .* maxTrans ;
        case 1
            clustCent(1:3,:) = clustCent(1:3,:) .* pi .*2;
            clustCent(4:6,:) = clustCent(4:6,:) .* boxSize;
        case 2
            clustCent(1:3,:) = clustCent(1:3,:) .* pi .* 2;
        case 3
            clustCent(1:3,:) = clustCent(1:3,:);
             clustCent(4,:) = clustCent(4,:);
        case 4
            clustCent(1:3,:) = clustCent(1:3,:) .* pi .* 2;
            %%Bug?
            clustCent(7,:) = clustCent(7,:) .* boxSize;
        case 5
            clustCent(1:3,:) = clustCent(1:3,:);
        case 6
            clustCent(1:3,:) = clustCent(1:3,:) .* pi .* 2;
    end
    
    clustCent = clustCent';
else
    if nargout > 2
        
        for cN = 1:numClust
            myMembers = find(data2cluster == cN);
            cluster2dataCell{cN} = myMembers;
        end
    end
end
toc


switch mode
    case 0
        save('D:\\dragon_clusterT7.mat','data2cluster')
    case 1
        save('D:\\dragon_clusterT6.mat','data2cluster')
    case 2
        save('D:\\baseballBat_clusterR.mat','data2cluster')
    case 3
        save('D:\\baseballBat_clusterP.mat','data2cluster')
    case 4
        save('D:\\baseballBat_clusterRP.mat','data2cluster')
    case 5
        save('D:\\baseballBat_clusterSV.mat','data2cluster')
    case 6
        save('D:\\baseballBat_clusterRV.mat','data2cluster')
end


