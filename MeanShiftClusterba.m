function [clustCent,data2cluster,cluster2dataMat,numClust,cluster2dataCell] = MeanShiftCluster(dataPts,mode,boxSize)




%*** Check input ****
% if nargin < 2
%     error('no bandwidth specified')
% end
% 
% if nargin < 3
%     plotFlag = true;
%     plotFlag = false;
% end
tic
if mode == 1
    d_max = max(dataPts',[],1);
    d_min = min(dataPts',[],1);
    L = d_max - d_min;
    bandWidth = max(L)/30;
else
    bandWidth = 0.15;
end
% save('C:\\Users\\WYZ95\\desktop\\test_cow.mat','dataPts')
plotFlag = false;

%**** Initialize stuff ***
[numDim,numPts] = size(dataPts);
numClust        = 0;
bandSq          = bandWidth^2;
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

while numInitPts

    tempInd         = ceil( (numInitPts-1e-6)*rand);        % pick a random seed point
    stInd           = initPtInds(tempInd);                  % use this point as start of mean
    myMean          = dataPts(:,stInd);                     % intilize mean to this points location
    myMembers       = [];                                   % points that will get added to this cluster                          
    thisClusterVotes    = zeros(1,numPts,'uint16');         % used to resolve conflicts on cluster membership

    while 1     %loop untill convergence
        metric = (repmat(myMean,1,numPts) - dataPts).^2;
%         sqDistToAll = sum(metric);    %dist squared from mean to all points still active
        if mode == 1
            sqDistToAll = Beta2 .* sum(metric(1:3,:)) + Beta3 .* sum(metric(4:6,:));    %dist squared from mean to all points still active
        else
            sqDistToAll = sum(metric);    %dist squared from mean to all points still active
        end
        inInds      = find(sqDistToAll < bandSq);               %points within bandWidth
        thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;  %add a vote for all the in points belonging to this cluster
        
        myOldMean   = myMean;                                   %save the old mean
        
        kernel_value = EpanechnikovKernel(numDim, size(inInds,2), sqDistToAll(inInds)./bandSq);
        numerator  = sum(dataPts(:,inInds) * kernel_value' , 2);
        denominator = sum(kernel_value , 2);
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
                pause
            end
        end

        %**** if mean doesn't move much stop this cluster ***
        if norm(myMean-myOldMean) < stopThresh
            
            %check for merge posibilities
            mergeWith = 0;
            for cN = 1:numClust
                distToOther = norm(myMean-clustCent(:,cN));     %distance from posible new clust max to old clust max
                if distToOther < bandWidth/2                    %if its within bandwidth/2 merge new and old
                    mergeWith = cN;
                    break;
                end
            end
            
            
            if mergeWith > 0    % something to merge
                clustCent(:,mergeWith)       = 0.5*(myMean+clustCent(:,mergeWith));             %record the max as the mean of the two merged (I know biased twoards new ones)
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

clusterSize = zeros(numClust,1);
[~,data2cluster] = max(clusterVotes,[],1);                %a point belongs to the cluster with the most votes
cluster2dataMat = -1 .* ones(numClust,numPts);
if mode == 1
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
    clustCent = clustCent';
else
    if nargout > 2
        cluster2dataCell = cell(numClust,1);
        for cN = 1:numClust
            myMembers = find(data2cluster == cN);
            cluster2dataCell{cN} = myMembers;
        end
    end
end
toc




