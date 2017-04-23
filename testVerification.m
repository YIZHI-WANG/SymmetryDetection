load vertex.mat
load face.mat

% v1 = 31949;
% v2 = 35925;
% transformation = [-9.85874 -13.7478 -5.00184];

v1 = 25345;
v2 = 30863;
transformation = [10.0232 15.1413 3.3433];

[ patch1, patch2 ] = Verification_entrance( v1, v2, transformation, vertex, face );
