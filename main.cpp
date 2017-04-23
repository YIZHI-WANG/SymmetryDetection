#ifndef IGL_VIEWER_WITH_NANOGUI
#define IGL_VIEWER_WITH_NANOGUI
#endif

#include "SymmetryDetector.h"

int main()
{
	////string file_path = SYMDETEC_SHARED_PATH "/cow.OFF";
	//double smooth_radius = 10;
	//double sample_ratio = 0.5;
	//double prune_threshold = 0.75;
	//double pair_threshold = 0.1;
	//int transformation_space_mode = TRANSFORM_6D;
	//SymmetryDetector symmetryDetector;
	//symmetryDetector.set_file_path(file_path);
	//symmetryDetector.set_smooth_radius(smooth_radius);
	//symmetryDetector.set_sample_ratio(sample_ratio);
	//symmetryDetector.set_prune_threshold(prune_threshold);
	//symmetryDetector.set_pair_threshold(pair_threshold);
	//symmetryDetector.set_transformation_space_mode(transformation_space_mode);
	//symmetryDetector.init_mesh();
	//symmetryDetector.init_viewer();
	//symmetryDetector.launch_viewer();
	
	SymmetryDetector symmetryDetector;
	symmetryDetector.init_viewer();
	symmetryDetector.launch_viewer();
}