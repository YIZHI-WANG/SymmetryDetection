#include "Mesh.h"
#include <igl/viewer/Viewer.h>
#include <igl/igl_inline.h>

class SymmetryDetector
{
public:
	SymmetryDetector();
	~SymmetryDetector();
	void set_file_path(string path);
	void set_remesh_fineness(double fineness);
	void set_smooth_radius(double radius);
	void set_sample_ratio(double ratio);
	void set_prune_threshold(double threshold);
	void set_pair_threshold(double threshold);
	void set_band_width(double bandwidth);
	void set_patching_threshold(double threshold);
	void set_transformation_space_mode(mode val);
	void init_mesh();
	void detect();
	void init_viewer();
	void launch_viewer();
private:
	Mesh* mesh;
	string file_path;
	double remesh_fineness;
	double smooth_radius;
	double pair_threshold;
	double prune_threshold;
	double sample_ratio;
	double band_width;
	double patching_threshold;
	bool do_patching;
	mode   transform_mode ;

	//the viewer body and parameter
	viewer::Viewer viewer;
	int current_show_cluster;
	int current_show_patch;
	int input_cluster_num;
	int input_patch_num;
	bool show_points;
	bool show_pair_lines;
	void set_viewer_parameter();
	void set_viewer_mesh();
	void set_viewer_menu();
	void set_viewer_keys();
	void clear_viewer();
	void clear_points();
	void clear_lines();
	void clear_points_and_lines();
	void show_sample_points();
	void show_pruned_points();
	void show_pairs();
	void show_cluster(int which_cluster);
	void show_patch(int which_cluster);
	void show_nth_patch(int n);
	void on_open_event();
};