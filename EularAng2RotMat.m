function [ RotMat ] = EularAng2RotMat( EularAng )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    RotMat = zeros(3,3);
	%%N
	RotMat(1, 1) = 1.0 * (cos(EularAng(1))*cos(EularAng(2)));
	RotMat(2, 1) = 1.0 * (sin(EularAng(1))*cos(EularAng(2)));
	RotMat(3, 1) = -1.0 * (sin(EularAng(2)));
	%%O
	RotMat(1, 2) = 1.0 * (cos(EularAng(1))*sin(EularAng(2))*sin(EularAng(3)) - sin(EularAng(1))*cos(EularAng(3)));
	RotMat(2, 2) = 1.0 * (sin(EularAng(1))*sin(EularAng(2))*sin(EularAng(3)) + cos(EularAng(1))*cos(EularAng(3)));
	RotMat(3, 2) = 1.0 * (cos(EularAng(2))*sin(EularAng(3)));
	%%A
	RotMat(1, 3) = 1.0 * (cos(EularAng(1))*sin(EularAng(2))*cos(EularAng(3)) + sin(EularAng(1))*sin(EularAng(3)));
	RotMat(2, 3) = 1.0 * (sin(EularAng(1))*sin(EularAng(2))*cos(EularAng(3)) - cos(EularAng(1))*sin(EularAng(3)));
	RotMat(3, 3) = 1.0 * (cos(EularAng(2))*cos(EularAng(3)));

end

