function [  ] = showResults(dataPts, data2cluster,mode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

metric = dataPts';

switch mode
    case 1
        [mappedX, mapping] = compute_mapping(metric, 'PCA', 2);
        x_max = max(mappedX(:,1));
        x_min = min(mappedX(:,1));
        y_max = max(mappedX(:,2));
        y_min = min(mappedX(:,2));
        clusterSize = max(data2cluster) + 1;
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';
        figure(1)
        plot(mappedX(:,1), mappedX(:,2),'o','markerfacecolor',cVec(mod(0,length(cVec)) + 1))
        axis([x_min x_max y_min  y_max]);
        title(['Before Clustering' ])
        figure(2)

        for i = 1:clusterSize
            thisCluster = find(data2cluster == i-1);
            plot(mappedX(thisCluster,1), mappedX(thisCluster,2),'o','markerfacecolor',cVec(mod(i,length(cVec)) + 1))
            axis([x_min x_max y_min  y_max]);

            hold on
    %         pause
    %         clf
        end
        title(['Clustering Result:' num2str(clusterSize)])
        figure(3)
        for i = 1:clusterSize
            thisCluster = find(data2cluster == i-1);
            plot(mappedX(thisCluster,1), mappedX(thisCluster,2),'o','markerfacecolor',cVec(mod(i,length(cVec)) + 1))
            axis([x_min x_max y_min  y_max]);
            title(['Cluster '  num2str(i-1)])
            hold on
            pause
            clf
        end
    case 2
        [mappedX, mappingx] = compute_mapping(metric(:,1:3), 'PCA', 1);
        [mappedY, mappingy] = compute_mapping(metric(:,4:6), 'PCA', 1);
        x_max = max(mappedX(:,1));
        x_min = min(mappedX(:,1));
        y_max = max(mappedY(:,1));
        y_min = min(mappedY(:,1));
        clusterSize = max(data2cluster) + 1;
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';
        figure(1)
        plot(mappedX(:,1), mappedY(:,1),'o','markerfacecolor',cVec(mod(0,length(cVec)) + 1))
        axis([x_min x_max y_min  y_max]);
        title(['Before Clustering' ])
        figure(2)

        for i = 1:clusterSize
            thisCluster = find(data2cluster == i-1);
            plot(mappedX(thisCluster,1), mappedY(thisCluster,1),'o','markerfacecolor',cVec(mod(i,length(cVec)) + 1))
            axis([x_min x_max y_min  y_max]);

            hold on
    %         pause
    %         clf
        end
        title(['Clustering Result:' num2str(clusterSize)])
        figure(3)
        for i = 1:clusterSize
            thisCluster = find(data2cluster == i-1);
            plot(mappedX(thisCluster,1), mappedY(thisCluster,1),'o','markerfacecolor',cVec(mod(i,length(cVec)) + 1))
            axis([x_min x_max y_min  y_max]);
            title(['Cluster '  num2str(i-1)])
            hold on
            pause
            clf
        end
    case 3
    figure
    scatter3(dataPts(1,:), dataPts(2,:), dataPts(3,:),3,data2cluster,'filled')
    title('Result of Cluster')
end
% figure
% scatter(mappedX(:,1), mappedX(:,2),3,data2cluster,'filled')
% title('Result of PCA')


end

