function [ vector_rot ] = Rot_By_Eular( vec, Eular )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
RotMat = EularAng2RotMat(Eular);
vector_rot = RotMat * vec;
vector_rot = vector_rot';
end

