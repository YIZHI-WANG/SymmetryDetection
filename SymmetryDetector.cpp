#include "SymmetryDetector.h"

#if !defined(__APPLE__)
std::string file_dialog(const std::vector<std::pair<std::string, std::string>> &filetypes, bool save) {
#define FILE_DIALOG_MAX_BUFFER 1024
#if defined(_WIN32)
	OPENFILENAME ofn;
	ZeroMemory(&ofn, sizeof(OPENFILENAME));
	ofn.lStructSize = sizeof(OPENFILENAME);
	char tmp[FILE_DIALOG_MAX_BUFFER];
	ofn.lpstrFile = tmp;
	ZeroMemory(tmp, FILE_DIALOG_MAX_BUFFER);
	ofn.nMaxFile = FILE_DIALOG_MAX_BUFFER;
	ofn.nFilterIndex = 1;

	std::string filter;

	if (!save && filetypes.size() > 1) {
		filter.append("Supported file types (");
		for (size_t i = 0; i < filetypes.size(); ++i) {
			filter.append("*.");
			filter.append(filetypes[i].first);
			if (i + 1 < filetypes.size())
				filter.append(";");
		}
		filter.append(")");
		filter.push_back('\0');
		for (size_t i = 0; i < filetypes.size(); ++i) {
			filter.append("*.");
			filter.append(filetypes[i].first);
			if (i + 1 < filetypes.size())
				filter.append(";");
		}
		filter.push_back('\0');
	}
	for (auto pair : filetypes) {
		filter.append(pair.second);
		filter.append(" (*.");
		filter.append(pair.first);
		filter.append(")");
		filter.push_back('\0');
		filter.append("*.");
		filter.append(pair.first);
		filter.push_back('\0');
	}
	filter.push_back('\0');
	ofn.lpstrFilter = filter.data();

	if (save) {
		ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;
		if (GetSaveFileNameA(&ofn) == FALSE)
			return "";
	}
	else {
		ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
		if (GetOpenFileNameA(&ofn) == FALSE)
			return "";
	}
	return std::string(ofn.lpstrFile);
#else
	char buffer[FILE_DIALOG_MAX_BUFFER];
	std::string cmd = "/usr/bin/zenity --file-selection ";
	if (save)
		cmd += "--save ";
	cmd += "--file-filter=\"";
	for (auto pair : filetypes)
		cmd += "\"*." + pair.first + "\" ";
	cmd += "\"";
	FILE *output = popen(cmd.c_str(), "r");
	if (output == nullptr)
		throw std::runtime_error("popen() failed -- could not launch zenity!");
	while (fgets(buffer, FILE_DIALOG_MAX_BUFFER, output) != NULL)
		;
	pclose(output);
	std::string result(buffer);
	result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
	return result;
#endif
}
#endif

SymmetryDetector::SymmetryDetector()
{
	mesh = NULL;
	//init some parameter
	set_remesh_fineness(5);
	set_smooth_radius(1000);
	set_sample_ratio(0.5);
	set_prune_threshold(0.85);
	set_pair_threshold(0.1);
	set_band_width(10);
	set_patching_threshold(0.05);
	set_transformation_space_mode(ROTATION_AND_SYMMETRY_7D);
	show_pair_lines = false;
}

SymmetryDetector::~SymmetryDetector()
{
}

void SymmetryDetector::set_file_path(string path)
{
	file_path = path;
}

void SymmetryDetector::set_remesh_fineness(double fineness)
{
	remesh_fineness = fineness;
}

void SymmetryDetector::set_smooth_radius(double radius)
{
	smooth_radius = radius;
}

void SymmetryDetector::set_sample_ratio(double ratio)
{
	sample_ratio = ratio;
}

void SymmetryDetector::set_prune_threshold(double threshold)
{
	prune_threshold = threshold;
}

void SymmetryDetector::set_pair_threshold(double threshold)
{
	pair_threshold = threshold;
}

void SymmetryDetector::set_band_width(double bandwidth)
{
	band_width = bandwidth;
}

void SymmetryDetector::set_patching_threshold(double threshold)
{
	patching_threshold = threshold;
}

void SymmetryDetector::set_transformation_space_mode(mode val)
{
	transform_mode = val;
}

void SymmetryDetector::init_mesh()
{
	mesh->readMesh(file_path);
	mesh->set_sample_ratio(sample_ratio);
	mesh->set_pair_threshold(pair_threshold);
	mesh->set_prune_threshold(prune_threshold);
	mesh->set_remesh_fineness(remesh_fineness);
	mesh->set_smooth_radius(smooth_radius);
	mesh->set_band_width(band_width);
	mesh->set_patching_threshold(patching_threshold);
	mesh->set_transformation_space_mode(transform_mode);
}

void SymmetryDetector::detect()
{
	clear_points_and_lines();
	set_viewer_parameter();
	cout << "Start Detecting..." << endl << endl;
	mesh->init_data();
	mesh->set_remesh_fineness(remesh_fineness);
	mesh->set_smooth_radius(smooth_radius);
	mesh->set_sample_ratio(sample_ratio);
	mesh->set_pair_threshold(pair_threshold);
	mesh->set_prune_threshold(prune_threshold);
	mesh->set_band_width(band_width);
	mesh->set_patching_threshold(patching_threshold);
	mesh->set_transformation_space_mode(transform_mode);
	mesh->remesh();
	mesh->compute_local_frame();
	mesh->sample();
	mesh->pairing();
	mesh->compute_transformation();
	mesh->meanShift_cluster();
	if (do_patching)
	{
		mesh->patching();
	}
	cout << "Finish Detecting!" << endl << endl;
}

void SymmetryDetector::init_viewer()
{
	set_viewer_parameter();
	do_patching = false;

	set_viewer_menu();

	set_viewer_keys();
}

void SymmetryDetector::launch_viewer()
{
	viewer.launch();
}

void SymmetryDetector::set_viewer_parameter()
{
	current_show_cluster = 0;
	current_patch = 0;
	current_show_patch = 0;
	input_cluster_num = 0;
	input_patch_num = 0;
	show_points = false;
}

void SymmetryDetector::clear_points()
{
	viewer.data.points.resize(0, 0);
}

void SymmetryDetector::clear_lines()
{
	viewer.data.lines.resize(0, 0);
}

void SymmetryDetector::clear_points_and_lines()
{
	clear_points();
	clear_lines();
}

void SymmetryDetector::set_viewer_mesh()
{
	viewer.data.set_mesh(mesh->get_vertex(), mesh->get_face());
	viewer.core.align_camera_center(viewer.data.V, viewer.data.F);
}

void SymmetryDetector::set_viewer_menu()
{
	// Extend viewer menu
	viewer.callback_init = [&](igl::viewer::Viewer& viewer)
	{
		viewer.ngui->setFixedSize(Eigen::Vector2i(120, 15));

		// Create nanogui widgets
		nanogui::Window *window = viewer.ngui->addWindow(Eigen::Vector2i(0, 0), "Symmetry Detector");

		
		viewer.ngui->addGroup("Detector Parameters");
		viewer.ngui->addVariable<double>("Sample Ratio", [&](double val) {
			sample_ratio = val; // set
		}, [&]() {
			return sample_ratio; // get
		});

		viewer.ngui->addVariable<double>("Remesh Fineness(mean_area/)", [&](double val) {
			remesh_fineness = val; // set
		}, [&]() {
			return remesh_fineness; // get
		});

		viewer.ngui->addVariable<double>("Smooth Radius(* boxSize)", [&](double val) {
			smooth_radius = val; // set
		}, [&]() {
			return smooth_radius; // get
		});

		viewer.ngui->addVariable<double>("Prune Threshold", [&](double val) {
			prune_threshold = val; // set
		}, [&]() {
			return prune_threshold; // get
		});

		viewer.ngui->addVariable<double>("Pair Threshold", [&](double val) {
			pair_threshold = val; // set
		}, [&]() {
			return pair_threshold; // get
		});

		viewer.ngui->addVariable<double>("Bandwidth", [&](double val) {
			band_width = val; // set
		}, [&]() {
			return band_width; // get
		});

		viewer.ngui->addVariable<double>("Patching Threshold", [&](double val) {
			patching_threshold = val; // set
		}, [&]() {
			return patching_threshold; // get
		});

		viewer.ngui->addVariable("Mode", transform_mode, true)
				   ->setItems({ "Transform_7D", "Transform_6D", "Rotation", "Reflection", "Rot & Ref", "Rot&Ref&Scal" /*, "Ref Vector", "Rot & RefVec"*/
							  });

		viewer.ngui->addVariable("Do Patching", do_patching);

		viewer.ngui->addGroup("Operations");
		// Add buttons
		viewer.ngui->addButton("Open", [&]() {
			on_open_event();
		});
		viewer.ngui->addButton("Detect", [&]() {
			detect();
		});

		viewer.ngui->addGroup("Show Options");
		viewer.ngui->addVariable<int>("Set Cluster", [&](int val) {
			input_cluster_num = val; // set
		}, [&]() {
			return input_cluster_num; // get
		});

		viewer.ngui->addVariable<int>("Set Patch", [&](int val) {
			input_patch_num = val; // set
		}, [&]() {
			return input_patch_num; // get
		});

		// Add buttons
		viewer.ngui->addButton("Show Sampling Points", [&]() {
			show_sample_points();
		});

		viewer.ngui->addButton("Show Pruned Points", [&]() {
			show_pruned_points();
		});

		viewer.ngui->addButton("Show Pairs", [&]() {
			show_pairs();
		});

		//viewer.ngui->addButton("Show/Close Samples", [&]() {
		//	show_points = !show_points;
		//});

		viewer.ngui->addButton("Show Set Cluster", [&]() {
			clear_points_and_lines();
			current_patch = 0;
			int cluster_num = mesh->get_cluster_num();

			if (input_cluster_num < 0)
			{
				cout << "Number should be >= 0." << endl;
			}
			else if (input_cluster_num >= cluster_num)
			{
				cout << "Number should be <= " << cluster_num - 1 << "." << endl;
			}
			else
			{
				show_cluster(input_cluster_num);
				current_show_cluster = input_cluster_num + 1;
				if (current_show_cluster >= cluster_num)
				{
					current_show_cluster = 0;
				}
			}
		});

		viewer.ngui->addButton("Show Next Cluster", [&]() {
			clear_points_and_lines();

			show_cluster(current_show_cluster);

			current_show_cluster++;
			current_patch = 0;
			if (current_show_cluster >= mesh->get_cluster_num())
			{
				current_show_cluster = 0;
			}
		});

		viewer.ngui->addButton("Show This Patch", [&]() {
			current_show_cluster--;
			if (current_show_cluster < 0)
			{
				cout << "I don't have this patch: " << current_show_cluster << "..." << endl;
			}
			else
			{
				int which_cluster = mesh->get_cluster_from_rank(current_show_cluster);
				int max_patch = mesh->get_max_patch(which_cluster);
				cout << "This is the ";
				switch (current_patch)
				{
				case 0:
					cout << "1st "; break;
				case 1:
					cout << "2nd "; break;
				case 2:
					cout << "3rd "; break;
				default:
					cout << current_patch + 1 << "th "; break;
				}
				cout << "patch of cluster " << current_show_cluster + 1 << "." << endl;
				show_patch(which_cluster, current_patch);
				current_patch++;
				if (current_patch == max_patch) current_patch = 0;
			}
			current_show_cluster++;
		});

		viewer.ngui->addButton("Show Set Patch", [&]() {
			clear_points_and_lines();

			int patch_num = mesh->get_patch_num();

			if (input_patch_num < 0)
			{
				cout << "Number should be >= 0." << endl;
			}
			else if (input_patch_num >= mesh->get_patch_num())
			{
				cout << "Number should be <= " << patch_num - 1 << "." << endl;
			}
			else
			{
				show_nth_patch(input_patch_num);
				current_show_patch = input_patch_num + 1;
				if (current_show_patch >= patch_num)
				{
					current_show_patch = 0;
				}
			}
		});

		viewer.ngui->addButton("Show Next Patch", [&]() {
			show_nth_patch(current_show_patch);

			current_show_patch++;
			
			if (current_show_patch >= mesh->get_patch_num())
			{
				current_show_patch = 0;
			}

		});

		viewer.ngui->addGroup("Draw Settings");
		viewer.ngui->addVariable("Show overlay", viewer.core.show_overlay);
		viewer.ngui->addVariable("Show overlay depth", viewer.core.show_overlay_depth);
		viewer.ngui->addVariable("Background", (nanogui::Color &) viewer.core.background_color);
		viewer.ngui->addVariable("Line color", (nanogui::Color &) viewer.core.line_color);
		viewer.ngui->addVariable("Wireframe", viewer.core.show_lines);
		viewer.ngui->addVariable("Fill", viewer.core.show_faces);
		viewer.ngui->addVariable("Show Pair Lines", show_pair_lines);
		viewer.ngui->addVariable("Show vertex labels", viewer.core.show_vertid);
		viewer.ngui->addVariable("Show faces labels", viewer.core.show_faceid);
		// Generate menu
		viewer.screen->performLayout();

		return false;
	};

}

void SymmetryDetector::set_viewer_keys()
{
	viewer.callback_key_down = [&](igl::viewer::Viewer& viewer, unsigned char key, int modifier)
	{
		//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
		if (key == 'S')
		{
			clear_points_and_lines();

			show_cluster(current_show_cluster);

			current_show_cluster++;
			if (current_show_cluster >= mesh->get_cluster_num())
			{
				current_show_cluster = 0;
			}
			current_patch = 0;
		}
		else if (key == 'P')
		{
			current_show_cluster--;
			if (current_show_cluster < 0)
			{
				cout << "I don't have this patch: " << current_show_cluster << "..." << endl;
			}
			else
			{
				int which_cluster = mesh->get_cluster_from_rank(current_show_cluster);
				int max_patch = mesh->get_max_patch(which_cluster);
				cout << "This is the ";
				switch (current_patch)
				{
				case 0:
					cout << "1st "; break;
				case 1:
					cout << "2nd "; break;
				case 2: 
					cout << "3rd "; break;
				default:
					cout << current_patch+1 << "th "; break;
				}
				cout << "patch of cluster " << current_show_cluster + 1 << "."<<endl;
				show_patch(which_cluster, current_patch);
				current_patch++;
				if (current_patch == max_patch) current_patch = 0;
			}
			current_show_cluster++;
			
			//current_show_cluster--;
			//if (current_show_cluster < 0)
			//{
			//	cout << "I don't have this patch: " << current_show_cluster << "..." << endl;
			//}
			//else
			//{
			//	//show_patch(current_show_cluster);
			//}
			//current_show_cluster++;
		}
		else if (key == 'N')
		{
			show_nth_patch(current_show_patch);

			current_show_patch++;
			if (current_show_patch >= mesh->get_patch_num())
			{
				current_show_patch = 0;
			}
		}
		return false;
	};
}

void SymmetryDetector::clear_viewer()
{
	viewer.data.clear();
}

void SymmetryDetector::show_sample_points()
{
	MatrixXd V = mesh->get_vertex();
	VectorXi V_Index_Sample = mesh->get_sample_index();
	viewer.data.points.resize(0, 0);

	for (int j = 0; j < V_Index_Sample.size(); j++)
	{
		viewer.data.add_points(V.row(V_Index_Sample(j)), Eigen::RowVector3d(0.25, 0.25, 1));
	}
}

void SymmetryDetector::show_pruned_points()
{
	MatrixXd V = mesh->get_vertex();
	VectorXi V_Index_Sample = mesh->get_pruned_index();
	viewer.data.points.resize(0, 0);

	for (int j = 0; j < V_Index_Sample.size(); j++)
	{
		viewer.data.add_points(V.row(V_Index_Sample(j)), Eigen::RowVector3d(0.25, 0.25, 1));
	}
}

void SymmetryDetector::show_pairs()
{
	MatrixXi Pair_Index = mesh->get_pairs();
	MatrixXd V = mesh->get_vertex();
	for (int i = 0; i < Pair_Index.rows(); i++)
	{
		viewer.data.add_edges
		(
			V.row(Pair_Index(i, 0)),
			V.row(Pair_Index(i, 1)),
			Eigen::RowVector3d(0.25, 0.25, 1)
		);
	}
}

void SymmetryDetector::show_cluster(int which_cluster)
{
	MatrixXd V = mesh->get_vertex();
	VectorXi vertex_in_cluster;
	which_cluster = mesh->get_cluster_from_rank(which_cluster);
	VectorXd clusterCenter = mesh->get_cluster_center(which_cluster);
	mesh->extract_vertex_from_cluster(which_cluster, vertex_in_cluster);
	int num = vertex_in_cluster.rows();
	MatrixXd V_in_Cluster;
	V_in_Cluster.resize(num, 3);
	for (int j = 0; j < num; j++)
	{
		V_in_Cluster.row(j) = V.row(vertex_in_cluster(j));
	}
	viewer.data.set_points(V_in_Cluster, Eigen::RowVector3d(1, 0, 0.5));
	for (int j = 0; j < num / 2; j++)
	{
		viewer.data.add_edges
		(
			V.row(vertex_in_cluster(j * 2)),
			V.row(vertex_in_cluster(j * 2 + 1)),
			Eigen::RowVector3d(0.25, 0.25, 1)
		);
	}
	cout << "Clsuter " << which_cluster << ": ";
	int dimension = clusterCenter.rows();
	for (int j = 0; j < dimension; j++)
	{
		cout << clusterCenter(j) << ", ";
	}
	cout << endl;
}

void SymmetryDetector::show_patch(int which_patch)
{
	//clear_points_and_lines();

	//show_cluster(current_show_cluster);
	//which_cluster = mesh->get_cluster_from_rank(which_cluster);
	MatrixXd V_in_patch, V;
	VectorXd Patch1, Patch2;
	V_in_patch.resize(1, 3);
	V = mesh->get_vertex();
	Patch1 = mesh->get_patch(1,which_patch);
	Patch2 = mesh->get_patch(2,which_patch);
	int v;
	int pointer;
	
	if (Patch1.rows() > 0)
	{
		pointer = 0;
		v = Patch1(pointer++);
		while (v >= 0)
		{
			V_in_patch = V.row(v);
			viewer.data.add_points(V_in_patch, Eigen::RowVector3d(0, 1, 0));
			v = Patch1(pointer++);
		}
	}
	
	if (Patch2.rows() > 0)
	{
		pointer = 0;
		v = Patch2(pointer++);
		while (v >= 0)
		{
			V_in_patch = V.row(v);
			viewer.data.add_points(V_in_patch, Eigen::RowVector3d(1, 0, 0));
			v = Patch2(pointer++);
		}
	}
	
}

void SymmetryDetector::show_patch(int which_cluster,int which_patch)
{
	//clear_points_and_lines();

	//show_cluster(current_show_cluster);
	//which_cluster = mesh->get_cluster_from_rank(which_cluster);
	MatrixXd V_in_patch, V;
	VectorXd Patch1, Patch2;
	V_in_patch.resize(1, 3);
	V = mesh->get_vertex();
	Patch1 = mesh->get_patch(1, which_cluster, which_patch);
	Patch2 = mesh->get_patch(2, which_cluster, which_patch);
	int v;
	int pointer;

	if (Patch1.rows() > 0)
	{
		pointer = 0;
		v = Patch1(pointer++);
		while (v >= 0)
		{
			V_in_patch = V.row(v);
			viewer.data.add_points(V_in_patch, Eigen::RowVector3d(0, 1, 0));
			v = Patch1(pointer++);
		}
	}

	if (Patch2.rows() > 0)
	{
		pointer = 0;
		v = Patch2(pointer++);
		while (v >= 0)
		{
			V_in_patch = V.row(v);
			viewer.data.add_points(V_in_patch, Eigen::RowVector3d(1, 0, 0));
			v = Patch2(pointer++);
		}
	}

}
void SymmetryDetector::show_nth_patch(int n)
{
	int which_patch = mesh->get_nth_patch_index(n);
	int which_cluster = mesh->get_cluster_from_patch(which_patch);
	int patch_size = mesh->get_patch_size(which_patch);
	VectorXd clusterCenter = mesh->get_cluster_center(which_cluster);
	MatrixXd V = mesh->get_vertex();
	int dimension = clusterCenter.rows();

	cout << "Patch " << current_show_patch << ", size = " << patch_size << "." << endl;
	cout << "Cluster " << which_cluster << ":";
	for (int j = 0; j < dimension; j++)
	{
		cout << clusterCenter(j) << ", ";
	}
	cout << endl;

	viewer.data.lines.resize(0, 0);
	viewer.data.points.resize(0, 0);
	VectorXi vertex_in_cluster;
	mesh->extract_vertex_from_cluster(which_cluster, vertex_in_cluster);
	int num = vertex_in_cluster.rows();
	MatrixXd V_in_Cluster;
	V_in_Cluster.resize(num, 3);
	for (int j = 0; j < num; j++)
	{
		V_in_Cluster.row(j) = V.row(vertex_in_cluster(j));
	}

	for (int j = 0; j < num / 2; j++)
	{
		VectorXd v1 = V.row(vertex_in_cluster(j * 2));
		VectorXd v2 = V.row(vertex_in_cluster(j * 2 + 1));
		stringstream l1, l2;
		l1 << vertex_in_cluster(j * 2);
		l2 << vertex_in_cluster(j * 2 + 1);
		if (show_pair_lines)
		{
			viewer.data.add_edges(v1.transpose(), v2.transpose(), Eigen::RowVector3d(0.25, 0.25, 1));
		}
		//viewer.data.add_label(v1.transpose(), l1.str());
		//viewer.data.add_label(v2.transpose(), l2.str());
	}

	show_patch(which_patch);

}

void SymmetryDetector::on_open_event()
{
	string new_file_path;
	new_file_path = file_dialog({ { "off", "Object File Format" },{ "obj","Alias Wavefront Object" } }, false);
	if (new_file_path != "")
	{
		file_path = new_file_path;
		if (mesh != NULL)
		{
			delete mesh;
		}
		clear_viewer();
		mesh = new Mesh();

		init_mesh();
		set_viewer_parameter();
		set_viewer_mesh();
	}
}

