#include "reference.h"

class Mesh
{
public:
	Mesh();
	~Mesh();
	void readMesh(string filePath);
	void set_remesh_fineness(double fineness);
	void set_smooth_radius(double radius);
	void set_sample_ratio(double ratio);
	void set_prune_threshold(double threshold);
	void set_pair_threshold(double threshold);
	void set_band_width(double bandwidth);
	void set_patching_threshold(double threshold);
	void set_transformation_space_mode(mode value);
	void init_data();
	void remesh();
	void compute_local_frame();
	void sample();
	void pairing();
	void compute_transformation();
	void meanShift_cluster();
	void extract_vertex_from_cluster(int which_cluster, VectorXi& V_in_Cluster);
	void patching();
	const MatrixXd& get_vertex();
	const MatrixXi& get_face();
	const MatrixXd& get_signature();
	const VectorXi& get_sample_index();
	const VectorXi& get_pruned_index();
	const MatrixXi& get_pairs();
	VectorXd get_cluster_center(int which_cluster);
	VectorXd get_patch(int which, int which_cluster, int which_patch);
	VectorXd get_patch(int which, int patch_row);
	int get_nth_patch_index(int n);
	int get_patch_size(int n);
	int get_cluster_num();
	int get_max_patch(int cluster);
	int get_cluster_from_rank(int rank);
	int get_cluster_from_patch(int which_patch);
	int get_patch_num();
	double get_boxSize();
private:
	MatrixXd vertex;
	MatrixXi face;

	int vertex_num, vertex_num_0;
	int face_num,face_num_0;

	double remesh_fineness;

	int sample_aim;

	double totalArea;
	double unit_area;

	double boxSize;
	double smooth_radius;
	MatrixXd pricinple_direction_min, 
		     pricinple_direction_max, 
		     local_normal;
	VectorXd pricinple_curvanture_min, 
			 pricinple_curvanture_max;
	VectorXd face_area;
	double min_area_0,max_area_0;

	MatrixXd signature;

	double sample_ratio;
	VectorXi sample_index, 
			 pruned_sample_index, 
			 sample_index_dictionary;
	int sample_num;

	double prune_threshold;
	
	VectorXd min_compare2max_sorted,
			 min_pricinple_val_sorted,
			 max_pricinple_val_sorted;
	VectorXi sample_index_sorted_min_compare2max,
			 sample_index_sorted_min_pricin_val,
			 sample_index_sorted_max_pricin_val;
	mode transformation_space_mode;
	double pair_threshold;
	MatrixXi pairs;
	int pair_num;

	MatrixXd transformation;

	double band_width;
	int cluster_num;
	MatrixXd cluster_center, data2cluster, cluster2data;
	MatrixXi cluster_size;
	VectorXi clusterRank;
	
	double patching_threshold;
	double patch_count;
	MatrixXd patch1, patch2, clsuter_2_patch;
	VectorXi patch_size, patch_rank;
	VectorXi patch_in_cluster_num, patch_2_clutser;

	// Matlab instance
	Engine* engine;

	void sort_vertex_in_face();
	void sort_vertex_in_face_by_row();
	void sort_vertex_in_face_by_col(int start, int end, int col);
	void compute_boxSize();
	void compute_curvature();
	void compute_signature();
	void compute_local_normal();
	void compute_area();
	void compute_min_area();
	void sample_from_facegroup(int samplePerUnit, int SampleFaceStart, 
							   int SampleFaceEnd, int BigFaceNum, 
							   const VectorXi& BigFaceIndex, int& sampleNum,
							   VectorXi& V_Sample);
	void prune_umbilic_points();
	void sort_sample_by_signature();
	void find_similar(int mode, const VectorXd& SV_col, const VectorXd& SV_Sorted, 
					 VectorXi& floor_index, VectorXi& celling_index);
	void start_engine();
	bool close_engine();
};

