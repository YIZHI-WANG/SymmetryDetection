function [ vector_rot ] = Rot_By_Eular( vec, Eular )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
RotMat = EularAng2RotMat(Eular);
vector_rot = RotMat * vec;
vector_rot = vector_rot';
end

