#ifndef IGL_VIEWER_WITH_NANOGUI
#define IGL_VIEWER_WITH_NANOGUI
#endif

#ifndef INF
#define INF 2147483647
#endif

#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <igl/sortrows.h>
#include <igl/find.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/matlab/MatlabWorkspace.h>
#include <igl/jet.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include "symDetec_shared_path.h"
using namespace igl;
using namespace Eigen;
using namespace std;

void RotateMatrix2EularAngle(const MatrixXd Rot, VectorXd& Eular);
void prune_umbilic_points(const MatrixXd& SV, const VectorXi& V_Index_Sample, VectorXi& V_Index_Prune, VectorXi& V_Index_Dic, double prune_threshold);
void compute_signature(const VectorXd& PV1, const VectorXd& PV2, MatrixXd& SV);
void samping(int size, double Sample_Rate, VectorXi& V_Index_Sample);
void sort_by_signature(const MatrixXd& SV, VectorXd& SV_Sorted_7d, VectorXd& SV_Sorted_PV1, VectorXd& SV_Sorted_PV2, VectorXi& V_Index_Sorted_7d, VectorXi& V_Index_Sorted_PV1, VectorXi& V_Index_Sorted_PV2);
void find_similar(const VectorXi& V_Index_Prune, const VectorXd& SV_col, const VectorXd& SV_Sorted, double pair_threshold, VectorXi& floor_index, VectorXi& celling_index);
void search_pair(int dimension, const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted_7d, const VectorXi& V_Index_Sorted_PV1, const VectorXi& V_Index_Sorted_PV2, const MatrixXd& SV, const VectorXd& SV_Sorted_7d, const VectorXd& SV_Sorted_PV1, const VectorXd& SV_Sorted_PV2, double pair_threshold, MatrixXi& Pair_Index);
void compute_local_frame(const MatrixXd& PD1, const MatrixXd& PD2, MatrixXd& LN);
void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int dimension, MatrixXd& Transformation);
int MeanShift_Cluster(const MatrixXd& dataPts, double bandwidth, MatrixXd& clusterCenter, MatrixXd& data2cluster, MatrixXd& cluster2data);
void Extract_Vertex(const MatrixXd& cluster2data, const MatrixXi& Pair_Index, int which_cluster, VectorXi& V_in_Cluster);

int main(int argc, char *argv[])
{
	double prune_threshold = 0.75;
	double pair_threshold = 0.05;
	int dimension = 6;
	double bandwidth = 0.75;

	MatrixXd V;
	MatrixXi F;

	MatrixXd PD1, PD2, LN;
	VectorXd PV1, PV2;
	MatrixXd SV;
	VectorXd SV_Sorted_PV1, SV_Sorted_PV2, SV_Sorted_7d;

	VectorXi V_Index_Sample, V_Index_Dic, V_Index_Prune, V_Index_Sorted_7d, V_Index_Sorted_PV1, V_Index_Sorted_PV2;
	MatrixXi Pair_Index;

	MatrixXd Transformation;

	int clusterNum;
	MatrixXd clusterCenter, data2cluster, cluster2data;

	// Load a mesh in OFF format
	readOBJ(SYMDETEC_SHARED_PATH "/plane (2).OBJ", V, F);
	cout << "READ OBJ FINISH!" << endl << endl;

	cout << "Model information:" << endl;
	cout << "Vertex: " << V.rows() << endl;
	cout << "Faces: " << F.rows() << endl;

	//Compute curvature directions via quadric fitting
	principal_curvature(V, F, PD1, PD2, PV1, PV2);
	cout << "COMPUTE PRINCIPAL CURVATURE FINISH!" << endl << endl;

	//Compute signature ¦Ò = PV1/PV2
	compute_signature(PV1, PV2, SV);
	cout << "COMPUTE SIGNATURE FINISH!" << endl << endl;
	
	compute_local_frame(PD1, PD2, LN);
	cout << "COMPUTE LOCAL FRAME FINISH!" << endl << endl;

	//Sampling
	samping(V.rows(), 0.5, V_Index_Sample);
	cout << "SAMPLE FINISH!" << endl;
	cout << "After Sampling, Vertex: " << V_Index_Sample.rows() << endl << endl;

	//Prune umbilic points
	prune_umbilic_points(SV, V_Index_Sample, V_Index_Prune, V_Index_Dic, prune_threshold);
	cout << "PRUNE FINISH!" << endl;
	cout << "After Pruning, Vertex: " << V_Index_Prune.rows() << endl << endl;

	//Sort Siganature
	sort_by_signature(SV, SV_Sorted_7d, SV_Sorted_PV1, SV_Sorted_PV2, V_Index_Sorted_7d, V_Index_Sorted_PV1, V_Index_Sorted_PV2);
	cout << "SORT FINISH!" << endl << endl;
	//for (int i = 0; i < SV.rows(); i++)
	//{
	//	cout << SV_Sorted_7d(i) <<" , " << SV(V_Index_Sorted_7d(i)) << endl;
	//}

	//Pairing
	search_pair(dimension, V_Index_Dic, V_Index_Prune, V_Index_Sorted_7d, V_Index_Sorted_PV1, V_Index_Sorted_PV2, SV, SV_Sorted_7d, SV_Sorted_PV1, SV_Sorted_PV2, pair_threshold, Pair_Index);
	cout << "PAIR FINISH!" << endl;
	cout << "After Pairing, Pair: " << Pair_Index.rows() << endl << endl;

	//Compute Transformation 
	compute_transformation(Pair_Index, V, PV1, PV2, PD1, PD2, LN, dimension, Transformation);
	cout << "COMPUTE TRANSFORMATION FINISH!" << endl << endl;

	//Mean-Shift
	clusterNum = MeanShift_Cluster(Transformation, bandwidth, clusterCenter, data2cluster, cluster2data);
	//data2cluster.transposeInPlace();
	cout << "MEAN-SHIFT FINISH!" << endl;
	cout << "After MEAN-SHIFT, cluster center number: " << clusterNum << endl;
	//for (int i = 0; i < clusterNum; i++)
	//{
	//	cout << "center" << i + 1 << ": ";
	//	for (int j = 0; j < clusterCenter.cols(); j++)
	//	{
	//		cout << clusterCenter(i, j) << ", ";
	//	}
	//	cout << endl;
	//}
	//cout << cluster2data.rows() << "," << cluster2data.cols() << endl;
	//for (int i = 0; i < data2cluster.rows(); i++)
	//{
	//	cout << "point " << i << " ---> " << data2cluster(i,0) << endl;
	//}
	cout << endl;

	//// Init the viewer
	viewer::Viewer viewer;

	VectorXi seed;
	MatrixXd C;
	seed.resize(clusterNum);
	for (int i = 0; i < clusterNum; i++)
	{
		seed(i) = i;
	}
	jet(seed, true, C);

	for (int i = 0; i < clusterNum; i++)
	{
		VectorXi vertex_in_cluster;
		Extract_Vertex(cluster2data, Pair_Index, i, vertex_in_cluster);
		int num = vertex_in_cluster.rows();
		MatrixXd V_in_Cluster;
		V_in_Cluster.resize(num, 3);
		for (int j = 0; j < num; j++)
		{
			V_in_Cluster.row(j) = V.row(vertex_in_cluster(j));
		}
		viewer.data.add_points(V_in_Cluster, C.row(i));
	}
	
	//// Plot the mesh
	viewer.data.set_mesh(V, F);

	viewer.launch();
	//
	//// Average edge length for sizing
	//const double avg = avg_edge_length(V, F);

	//// Draw a blue segment parallel to the minimal curvature direction
	//const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
	//viewer.data.add_edges(V + PD1*avg, V - PD1*avg, blue);

	//// Draw a red segment parallel to the maximal curvature direction
	//viewer.data.add_edges(V + PD2*avg, V - PD2*avg, red);

	//// Hide wireframe
	//viewer.core.show_lines = false;
	//

}

void prune_umbilic_points(const MatrixXd& V, const MatrixXd& PD1, const MatrixXd& PD2, const VectorXd& PV1, const VectorXd& PV2, VectorXi& V_Index_Prune, double prune_threshold)
{
	VectorXi _V_Index_Prune;
	int size = V.rows();
	int size_P = 0;
	_V_Index_Prune.resize(size, 1);
	
	for (int i = 0; i < size; i++) 
	{
		double _pv1 = PV1(i);
		double _pv2 = PV2(i);
		//cout << i << "," << "k1 = " << _pv1 << ", k2 = " << _pv2 << endl;
		if (_pv2 == 0.0 && _pv1 != 0.0)
		{
			_V_Index_Prune(size_P++) = i;
			//cout << "0 add " << i << endl;
		}
		else if (fabs(_pv1 / _pv2) < prune_threshold)
		{
			_V_Index_Prune(size_P++) = i;
			//cout << "add " << i << endl;
		}
	}
	V_Index_Prune = _V_Index_Prune.head(size_P);
}

void prune_umbilic_points(const MatrixXd& SV, const VectorXi& V_Index_Sample, VectorXi& V_Index_Prune, VectorXi& V_Index_Dic, double prune_threshold)
{
	VectorXi _V_Index_Prune;
	int size = SV.rows();
	int size_Sample = V_Index_Sample.rows();
	int size_P = 0;
	_V_Index_Prune.resize(size_Sample, 1);
	V_Index_Dic.resize(size, 1);
	V_Index_Dic.fill(-1);

	for (int i = 0; i < size_Sample; i++)
	{
		//cout << i << "SV = " << SV(i) << endl;
		if (fabs(SV(V_Index_Sample(i),0)) < prune_threshold)
		{
			V_Index_Dic(V_Index_Sample(i)) = size_P;
			_V_Index_Prune(size_P++) = V_Index_Sample(i);
			//cout << "add " << i << endl;
		}
	}

	V_Index_Prune = _V_Index_Prune.head(size_P);
}

void compute_signature(const VectorXd& PV1, const VectorXd& PV2, MatrixXd& SV)
{
	int size = PV1.rows();
	SV.resize(size, 3);

	for (int i = 0; i < size; i++)
	{
		if(PV2(i) == 0.0)
		{
			SV(i,0) = (PV1(i) > 0.0) ? (-INF) : (INF);
		}
		else
		{
			SV(i,0) = PV1(i) / PV2(i);
		}
		//cout << "SV" << i << " = " << SV(i) << endl;
	}

	SV.col(1) = PV1;
	SV.col(2) = PV2;

}

void samping(int size, double Sample_Rate, VectorXi& V_Index_Sample)
{
	int sample_size = size * Sample_Rate;
	V_Index_Sample.resize(sample_size);

	vector<int> s_stl;
	for (int i = 0; i<size; ++i) s_stl.push_back(i);
	random_shuffle(s_stl.begin(), s_stl.end());

	for (int j = 0; j < sample_size; j++)
	{
		V_Index_Sample(j) = s_stl.at(j);
	}

}

void sort_by_signature(const MatrixXd& SV, VectorXd& SV_Sorted_7d, VectorXd& SV_Sorted_PV1, VectorXd& SV_Sorted_PV2, VectorXi& V_Index_Sorted_7d, VectorXi& V_Index_Sorted_PV1, VectorXi& V_Index_Sorted_PV2)
{
	VectorXd SV_7d = SV.col(0);
	VectorXd SV_PV1 = SV.col(1);
	VectorXd SV_PV2 = SV.col(2);
	sortrows(SV_7d, true, SV_Sorted_7d, V_Index_Sorted_7d);
	sortrows(SV_PV1, true, SV_Sorted_PV1, V_Index_Sorted_PV1);
	sortrows(SV_PV2, true, SV_Sorted_PV2, V_Index_Sorted_PV2);
}

void find_similar(const VectorXi& V_Index_Prune, const VectorXd& SV_col, const VectorXd& SV_Sorted, double pair_threshold, VectorXi& floor_index, VectorXi& celling_index)
{
	int pruned_point_size = V_Index_Prune.rows();
	int point_size = SV_col.rows();
	floor_index.resizeLike(V_Index_Prune);
	celling_index.resizeLike(V_Index_Prune);
	for (int i = 0; i < pruned_point_size; i++)
	{
		floor_index(i) = -1;
		celling_index(i) = -1;
		double _sv = SV_col(V_Index_Prune(i));
		double _sv_floor, _sv_celling;
		if (_sv > 0)
		{
			_sv_floor = _sv * (1 - pair_threshold);
			_sv_celling = _sv * (1 + pair_threshold);
		}
		else
		{
			_sv_floor = _sv * (1 + pair_threshold);
			_sv_celling = _sv * (1 - pair_threshold);
		}

		for (int p = 0; p < point_size; p++) {
			double _sv_current = SV_Sorted(p);
			if (_sv_current >= _sv_floor)
			{
				floor_index(i) = p;
				break;
			}
		}
		for (int q = point_size - 1; q >= 0; q--) {
			double _sv_current = SV_Sorted(q);
			if (_sv_current <= _sv_celling)
			{
				celling_index(i) = q;
				break;
			}
		}
	}
}

void search_pair(int dimension, const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted_7d, const VectorXi& V_Index_Sorted_PV1, const VectorXi& V_Index_Sorted_PV2, const MatrixXd& SV, const VectorXd& SV_Sorted_7d, const VectorXd& SV_Sorted_PV1, const VectorXd& SV_Sorted_PV2, double pair_threshold, MatrixXi& Pair_Index)
{
	int point_size = SV.rows();
	int pruned_point_size = V_Index_Prune.rows();
	int max_size = pruned_point_size * (pruned_point_size - 1) / 2;
	int pair_size = 0;

	MatrixXi _Pair_Index(max_size, 2);

	if (dimension == 7)
	{
		VectorXi floor_index, celling_index;
		find_similar(V_Index_Prune, SV.col(0), SV_Sorted_7d, pair_threshold, floor_index, celling_index);
		for (int i = 0; i < pruned_point_size; i++)
		{
			/*double _sv = SV(V_Index_Prune(i), 0);
			double _sv_floor, _sv_celling;
			if (_sv > 0)
			{
				_sv_floor = _sv * (1 - pair_threshold);
				_sv_celling = _sv * (1 + pair_threshold);
			}
			else
			{
				_sv_floor = _sv * (1 + pair_threshold);
				_sv_celling = _sv * (1 - pair_threshold);
			}

			int floor_index = -1, celling_index = -1;
			for (int p = 0; p < point_size; p++) {
				double _sv_current = SV_Sorted_7d(p);
				if (_sv_current >= _sv_floor)
				{
					floor_index = p;
					break;
				}
			}
			for (int q = point_size - 1; q >= 0; q--) {
				double _sv_current = SV_Sorted_7d(q);
				if (_sv_current <= _sv_celling)
				{
					celling_index = q;
					break;
				}
			}*/
			for (int k = floor_index(i); k <= celling_index(i); k++)
			{
				if (V_Index_Dic(V_Index_Sorted_7d(k)) != -1 && V_Index_Sorted_7d(k) != V_Index_Prune(i))
				{
					_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted_7d(k);
				}
			}
		}
	}
	else if (dimension == 6)
	{
		VectorXi floor_index_PV1, floor_index_PV2, celling_index_PV1, celling_index_PV2;
		find_similar(V_Index_Prune, SV.col(1), SV_Sorted_PV1, pair_threshold, floor_index_PV1, celling_index_PV1);
		find_similar(V_Index_Prune, SV.col(2), SV_Sorted_PV2, pair_threshold, floor_index_PV2, celling_index_PV2);
		for (int i = 0; i < pruned_point_size; i++)
		{	
			for (int m = floor_index_PV1(i); m <= celling_index_PV1(i); m++)
			{
				int V_Index_Sorted_PV1_ = V_Index_Sorted_PV1(m);
				for (int n = floor_index_PV2(i); n <= celling_index_PV2(i); n++)
				{
					int V_Index_Sorted_PV2_ = V_Index_Sorted_PV2(n);
					if (V_Index_Sorted_PV1_ == V_Index_Sorted_PV2_)
					{
						if (V_Index_Dic(V_Index_Sorted_PV1_) != -1 && V_Index_Sorted_PV1_ != V_Index_Prune(i))
						{
							//cout << SV(V_Index_Prune(i), 1) << "," << SV(V_Index_Prune(i), 2) << " and "<<SV(V_Index_Sorted_PV1_, 1) << "," << SV(V_Index_Sorted_PV1_, 2) << endl;
							_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted_PV1_;
						}
						break;
					}
				}
			}
		}
	}
	
	//cout << "pair_size = " << pair_size << endl;
	Pair_Index = _Pair_Index.block(0, 0, pair_size, 2);
}

void compute_local_frame(const MatrixXd& PD1, const MatrixXd& PD2, MatrixXd& LN)
{
	int size = PD1.rows();
	LN.resize(size, 3);
	for (int i = 0; i < size; i++)
	{
		Vector3d d1 = PD1.row(i);
		Vector3d d2 = PD2.row(i);
		LN.row(i) = d1.cross(d2);
	}
}

void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int dimension, MatrixXd& Transformation)
{
	int size = Pair_Index.rows();
	Transformation.resize(size, dimension);
	
	for (int i = 0; i < size; i++)
	{
		MatrixXd Frame_p1, Frame_p2;
		MatrixXd Rot;
		VectorXd Eular;
		VectorXd Trans;
		double Scale;
		Frame_p1.resize(3, 3);
		Frame_p2.resize(3, 3);
		Rot.resize(3, 3);
		Trans.resize(3);

		int pair_p1 = Pair_Index(i, 0);
		int pair_p2 = Pair_Index(i, 1);

		if (dimension == 7)
		{
			Scale = (PV1(pair_p1) / PV1(pair_p2) + PV2(pair_p1) / PV2(pair_p2)) / 2;
		}
		else if (dimension == 6)
		{
			Scale = 1;
		}

		Frame_p1.col(0) = LN.row(pair_p1).transpose();
		Frame_p1.col(1) = PD1.row(pair_p1).transpose();
		Frame_p1.col(2) = PD2.row(pair_p1).transpose();

		Frame_p2.col(0) = LN.row(pair_p2).transpose();
		Frame_p2.col(1) = PD1.row(pair_p2).transpose();
		Frame_p2.col(2) = PD2.row(pair_p2).transpose();

		Rot = Frame_p1 * Frame_p2.inverse();
		RotateMatrix2EularAngle(Rot, Eular);
		
		Trans = V.row(pair_p2).transpose() - Scale * Rot * (V.row(pair_p1).transpose());

		//cout << "point1 = " << V(pair_p1, 0) << "," << V(pair_p1, 1) << "," << V(pair_p1, 2) << endl;
		//cout << "point2 = " << V(pair_p2, 0) << "," << V(pair_p2, 1) << "," << V(pair_p2, 2) << endl;

		//cout << "Scale = " << Scale << endl;

		//cout << "Rot = " << endl
		//	<< Rot(0, 0) << "," << Rot(0, 1) << "," << Rot(0, 2) << "," << endl
		//	<< Rot(1, 0) << "," << Rot(1, 1) << "," << Rot(1, 2) << "," << endl
		//	<< Rot(2, 0) << "," << Rot(2, 1) << "," << Rot(2, 2) << "," << endl;

		//cout << "Eular = " << Eular(0) << "," << Eular(1) << "," << Eular(2) << endl;

		//cout << "Trans = " << Trans(0) << "," << Trans(1) << "," << Trans(2) << endl;

		//cout << endl;

		if (dimension == 7)
		{
			Transformation.row(i) << Scale, Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
		}
		else if (dimension == 6)
		{
			Transformation.row(i) << Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
		}
	}
}

void RotateMatrix2EularAngle(const MatrixXd Rot, VectorXd& Eular)
{
	Eular.resize(3);
	double 
	r1 = Rot(0, 0),  r2 = Rot(0, 1),   r3 = Rot(0, 2),
	r4 = Rot(1, 0),  r5 = Rot(1, 1),   r6 = Rot(1, 2),
	r7 = Rot(2, 0),  r8 = Rot(2, 1),   r9 = Rot(2, 2);
	double Rz = atan(r4 / r1);
	double sz = sin(Rz), cz = cos(Rz);
	double Ry = atan((-r7) / (r1*cz + r4*sz));
	double Rx = atan((r3*sz - r6*cz) / (r5*cz - r2*sz));

	Eular << Rx, Ry, Rz;
}

int MeanShift_Cluster(const MatrixXd& dataPts, double bandwidth, MatrixXd& clusterCenter, MatrixXd& data2cluster, MatrixXd& cluster2data)
{
	// Matlab instance
	Engine* engine;

	// Launch MATLAB
	matlab::mlinit(&engine);

	// Send Transformation matrix to matlab
	MatrixXd dataPts_ = dataPts.transpose();
	matlab::mlsetmatrix(&engine, "dataPts", dataPts_);

	// Send bandwidth to matlab
	//matlab::mlsetscalar(&engine, "bandwidth", bandwidth);

	//Cluster using matlab
	igl::matlab::mleval(&engine, "[clustCent,data2cluster,cluster2data,numClust] = MeanShiftCluster(dataPts)");
	
	// Get clustering
	matlab::mlgetmatrix(&engine, "clustCent", clusterCenter);
	matlab::mlgetmatrix(&engine, "data2cluster", data2cluster);
	matlab::mlgetmatrix(&engine, "cluster2data", cluster2data);
	int numClust = (int)matlab::mlgetscalar(&engine, "numClust");

	//cout << "matlab" << matlab::mleval(&engine, "size(clustCent)") << endl;
	//cout << "matlab" << matlab::mleval(&engine, "size(data2cluster)") << endl;
	//cout << "matlab" << matlab::mleval(&engine, "cluster2data(1,:)") << endl;
	//cout << data2cluster.rows() << " , " << data2cluster.cols() << endl;
	//cout << "Recieve:" << endl;
	//for (int i = 0; i < cluster2data.cols(); i++)
	//{
	//	cout << cluster2data(0, i) << " ";
	//}
	//cout << endl;
	//cout << cluster2data.rows() << " , " << cluster2data.cols() << endl;

	return numClust;
}

void Extract_Vertex(const MatrixXd& cluster2data, const MatrixXi& Pair_Index, int which_cluster, VectorXi& V_in_Cluster)
{
	VectorXi Pair_cluster;
	Pair_cluster.resize(cluster2data.cols());
	int num = 0;
	for (int i = 0; i < cluster2data.cols(); i++)
	{
		if (cluster2data(which_cluster, i) >= 0)
		{
			Pair_cluster(i) = cluster2data(which_cluster, i);
			num++;
		}
	}

	V_in_Cluster.resize(num * 2);
	for (int i = 0; i < num; i++)
	{
		//TODO CHECK REPEAT VERTEX
		V_in_Cluster(i * 2) = Pair_Index(Pair_cluster(i), 0);
		V_in_Cluster(i * 2 + 1) = Pair_Index(Pair_cluster(i), 1);
	}
}