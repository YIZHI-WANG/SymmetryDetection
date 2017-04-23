#include <igl/viewer/Viewer.h>
using namespace igl;

class M_Viewer
{
public:
	M_Viewer();
	~M_Viewer();
	void init_viewer();

private:
	//the viewer body and parameter
	viewer::Viewer viewer;
	int input_cluster;

};

