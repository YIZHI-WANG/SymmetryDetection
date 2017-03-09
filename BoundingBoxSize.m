function [ boxSize ] = BoundingBoxSize( dataPts )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
maxPos          = max(dataPts,[],2);                          %biggest size in each dimension
minPos          = min(dataPts,[],2);                          %smallest size in each dimension
boundBox        = maxPos-minPos;                        %bounding box size
boxSize       = norm(boundBox);                       %indicator of size of data space

end

