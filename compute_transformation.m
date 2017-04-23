function [ scale, rot, trans ] = compute_transformation( mode, n1,n2,pd11,pd12,pd21,pd22,pv11,pv12,pv21,pv22,pos1,pos2 )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

    MatrixXd Frame_p1, Frame_p2;
    rot1 = zeros(3,3);
    rot2 = zeros(3,3);


    if mode == 0
		scale = (pv11 / pv12+ pv21 / pv22) / 2;
    else
        scale = 1;
    end
    
	origin = [0 0 0];

    p1_p2_n = n1' * n2 ;
    p1_p2_n = p1_p2_n ./ norm(p1_p2_n);
	cos_n1_n2 = n1 * n2';
	angle_n1_n2 = acos(cos_n1_n2);
	rot1 = compute_rotation(origin, p1_p2_n, angle_n1_n2);
    
	pd11_rot1 = rot1' * pd11;

    cos_p1_p2 = pd11 * pd12' ;
    if cos_p1_p2 > 1.0
        cos_p1_p2 = 1.0;
    end
    angle_p1_p2 = acos(cos_p1_p2);
	normal = pd11_rot1' * pd12;

    if normal*n1' < 0
			angle_p1_p2 = -angle_p1_p2;
    end
    
    rot2 = compute_rotation(origin, N1_ROT, angle_p1_p2);
    
	rot = rot2 * rot1;
    
    trans = pos2 - po1;

end

