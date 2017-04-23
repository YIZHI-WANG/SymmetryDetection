function [ dist ] = compute_dist_of_point2tri( vertex, face, point, tri )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
v0 = vertex(:,face(1,tri));
v1 = vertex(:,face(2,tri));
v2 = vertex(:,face(3,tri));
% point = vertex(:,point);

B = v0;
E0 = v1 - B;
E1 = v2 - B;
D = B - point;

a = dot(E0,E0);
b = dot(E0,E1);
c = dot(E1,E1);
d = dot(E0,D);
e = dot(E1,D);
f = dot(D,D);

det = a*c - b*b;
s   = b*e - c*d;
t   = b*d - a*e;

if (s+t) <= det
	if s < 0
        %%region 4
		if t < 0 
            if d < 0
                t = 0;
                if (-d) >= a
                    s = 1;
					sqrDistance = a + 2*d +f;
                else
                    s = -d / a ;
					sqrDistance = d*s + f;
                end
            else
                s = 0;
                if e >= 0
                    t = 0;
                    sqrDistance = f;
                else
					if (-e) >= c
                        t = 1;
                        sqrDistance = c + 2*e +f;
                    else
                        t = -e / c;
                        sqrDistance = e*t + f;
                    end
                end
            end
        else %% region 3
			s = 0;
            if e >= 0
                t = 0;
                sqrDistance = f;
            else
				if ((-e) > c)
                    t = 1;
                    sqrDistance = c + 2*e + f;
                else
					t = -e / c ;
                    sqrDistance = e*t + f;
                end  
            end
        end
    else      
        if (t < 0)
            %% region 5
            t = 0;
            if ( d >= 0)
                s = 0;
                sqrDistance = f;
            else
				if ((-d) > a)
					s = 1;
                    sqrDistance = a + d*s + f;
                else
					s = -d / a;
                    sqrDistance = d*s + f;
                end
			end%% of region 5
        else
            %% region 0 
            invDet = 1 / det;
            s = s * invDet;
            t = t * invDet;
            sqrDistance = s * (a*s + b*t + 2*d) + t * (b*s + c*t + 2*e) + f;
        end
    end
else
	if (s < 0)
		 %% region 2
        tmp0 = b + d;
        tmp1 = c + e;
        if (tmp1 > tmp0)
			numer = tmp1 - tmp0;
			denom = a - 2*b + c;
            if (numer >= denom)
				s = 1;
				t = 0;
				sqrDistance = a + 2*b +f;
			else
				s = numer / denom;
                t = 1 - s;
                sqrDistance = s * (a*s + b*t + 2*d)+ t * (b*s + c*t + 2*e) + f;
            end
        else
			s = 0;
			if (tmp1 <= 0)
				t = 1;
                sqrDistance = c + 2*e + f;
            else
				if (e >= 0)
					t = 0;
					sqrDistance = f;
				else
					t = -e / c;
					sqrDistance = e*t + f;
                end
            end
        end %% of region 2
    else
		if (t < 0)
            %% region 6
            tmp0 = b + e;
            tmp1 = a + d;
            if (tmp1 > tmp0)
				numer = tmp1  - tmp0;
                denom = a - 2*b + c;
                if (numer >= denom)
					t = 1;
                    s = 0;
                    sqrDistance = c + 2*e + f;
                else
                    t = numer / denom;
                    s = 1 - t;
                    sqrDistance = s * (a*s + b*t + 2*d)+ t * (b*s + c*t + 2*e) + f;
                end
            else
				t = 0;
                if (tmp1 < 0)
					s = 1;
                    sqrDistance = a + 2*d + f;
                else
					if (d >= 0)
						s = 0;
                        sqrDistance = f;
                    else
						s = -d / a;
                        sqrDistance = d*s + f;
                    end
                end
			end%% of region 6
        else
			%% region 1
            numer = c + e - b - d;
            if (numer <= 0)
				s = 0;
                t = 1;
                sqrDistance = c + 2*e + f;
            else
				denom = a - 2*b + c;
                if (numer >= denom)
					s = 1;
                    t = 0;
                    sqrDistance = a + 2*d + f;
                else
					s = numer / denom;
                    t = 1 - s;
                    sqrDistance = s * (a*s + b*t + 2*d)+ t * (b*s + c*t + 2*e) + f;
                end
            end%% of region 1
        end
    end
end


if (sqrDistance < 0)
    sqrDistance = 0;
end

dist = sqrt(sqrDistance);
    
end

