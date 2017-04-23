#include <Eigen/Geometry>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/sortrows.h>
#include <igl/doublearea.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/matlab/MatlabWorkspace.h>
#include <algorithm>
#include <iostream>
using namespace Eigen;
using namespace std;
using namespace igl;

#ifndef INF
#define INF 2147483647
#endif

#ifndef SYMDETEC_SHARED_PATH
#define SYMDETEC_SHARED_PATH "../shared"
#endif


//#define TRANSFORM_7D		     0
//#define TRANSFORM_6D	         1
//#define ROTATION_3D	             2
//#define SYMMETRY_PLANE_4D        3
//#define ROTATION_AND_SYMMETRY_7D 4
//#define SYMMETRY_VECTOR_3D       5
//#define ROTATION_AND_SVECTOR_6D  6

enum mode {
	TRANSFORM_7D = 0,
	TRANSFORM_6D = 1,
	ROTATION_3D = 2,
	SYMMETRY_PLANE_4D = 3,
	ROTATION_AND_SYMMETRY_7D = 4,
	ROT_AND_SYM_AND_SCAL = 5,
	ROTATION_AND_SVECTOR_6D = 6
};

