#include "M_Viewer.h"

int temp;

M_Viewer::M_Viewer()
{
}

M_Viewer::~M_Viewer()
{
}

void M_Viewer::init_viewer()
{
	// Extend viewer menu
	viewer.callback_init = [&](igl::viewer::Viewer& viewer)
	{
		viewer.ngui->addVariable<int>("Show Which Cluster", [&temp](int val) {
			temp = val; // set
		}, [&]() {
			return input_cluster; // get
		});
		// Add buttons
		viewer.ngui->addButton("Show Sampling Points", [V, V_Index_Sample]() {
			viewer.data.points.resize(0, 0);

			for (int j = 0; j < V_Index_Sample.size(); j++)
			{
				viewer.data.add_points(V.row(V_Index_Sample(j)), Eigen::RowVector3d(0, 0, 1));
			}
		});
	};
}
