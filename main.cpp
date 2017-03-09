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
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include "symDetec_shared_path.h"
using namespace igl;
using namespace Eigen;
using namespace std;

void RotateMatrix2EularAngle(const MatrixXd& Rot, VectorXd& Eular);
void RotateEularAngle2Matrix(const VectorXd& Eular, MatrixXd& Rot);
void prune_umbilic_points(const MatrixXd& SV, const VectorXi& V_Index_Sample, VectorXi& V_Index_Prune, VectorXi& V_Index_Dic, double prune_threshold);
void compute_signature(const VectorXd& PV1, const VectorXd& PV2, MatrixXd& SV);
void samping(int size, double Sample_Rate, VectorXi& V_Index_Sample);
void sort_by_signature(const MatrixXd& SV, VectorXd& SV_Sorted_7d, VectorXd& SV_Sorted_PV1, VectorXd& SV_Sorted_PV2, VectorXi& V_Index_Sorted_7d, VectorXi& V_Index_Sorted_PV1, VectorXi& V_Index_Sorted_PV2);
void find_similar(const VectorXi& V_Index_Prune, const VectorXd& SV_col, const VectorXd& SV_Sorted, double pair_threshold, VectorXi& floor_index, VectorXi& celling_index);
void search_pair(int dimension, const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted_7d, const VectorXi& V_Index_Sorted_PV1, const VectorXi& V_Index_Sorted_PV2, const MatrixXd& SV, const VectorXd& SV_Sorted_7d, const VectorXd& SV_Sorted_PV1, const VectorXd& SV_Sorted_PV2, double pair_threshold, MatrixXi& Pair_Index);
void compute_local_frame(const MatrixXd& PD1, const MatrixXd& PD2, MatrixXd& LN);
void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int dimension, MatrixXd& Transformation);
void compute_rotation(const VectorXd& V_AXIS_start, const VectorXd& V_AXIS_end, const double angle, MatrixXd& Rot, MatrixXd& quaternionRot);
int MeanShift_Cluster(const MatrixXd& V, const MatrixXd& dataPts, MatrixXd& clusterCenter, MatrixXd& data2cluster, MatrixXd& cluster2data);
void Extract_Vertex(const MatrixXd& cluster2data, const MatrixXi& Pair_Index, int which_cluster, VectorXi& V_in_Cluster);
void Show_Pair_Points(const MatrixXd& V, const MatrixXi& Pair_Index, int num);
void printVector(const VectorXi& vec);
void printVector(const VectorXd& vec);
void printMatrix(const MatrixXd& mat);
void printMatrix(const MatrixXi& mat);

//// Init the viewer
viewer::Viewer m_viewer;

int main(int argc, char *argv[])
{
	double prune_threshold = 0.75;
	double pair_threshold = 0.05;
	int dimension = 6;

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
	//readOBJ(SYMDETEC_SHARED_PATH "/bunny.OBJ", V, F);
	//readOFF(SYMDETEC_SHARED_PATH "/cow.OFF", V, F);
	//readOFF(SYMDETEC_SHARED_PATH "/cheburashka.off", V, F);
	//readOBJ(SYMDETEC_SHARED_PATH "/cheburashka_arm.OBJ", V, F);
	//readOFF(SYMDETEC_SHARED_PATH "/decimated-knight.off", V, F);
	readOBJ(SYMDETEC_SHARED_PATH "/dummy.OBJ", V, F);
	//readOBJ(SYMDETEC_SHARED_PATH "/chair_005.OBJ", V, F);
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
	samping(V.rows(), 0.1, V_Index_Sample);
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
	clusterNum = MeanShift_Cluster(V,Transformation, clusterCenter, data2cluster, cluster2data);
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
	
	//// Plot the mesh
	m_viewer.data.set_mesh(V, F);

	int cn = 0;
	int pair_n = 0;
	// Extend viewer menu
	m_viewer.callback_init = [&](igl::viewer::Viewer& viewer)
	{
		// Add a button
		m_viewer.ngui->addButton("Show Pair", [V, Pair_Index, &pair_n]() {
			/*for (int i = 0; i < Pair_Index.rows(); i++)
			{
				m_viewer.data.add_edges
				(
					V.row(Pair_Index(i, 0)),
					V.row(Pair_Index(i, 1)),
					Eigen::RowVector3d(1, 0, 0)
				);
			}*/
			Show_Pair_Points(V, Pair_Index, -1);
			pair_n++;
			if (pair_n >= Pair_Index.rows())
			{
				pair_n = 0;
			}
		});

		m_viewer.ngui->addButton("Show Cluster", [V, clusterCenter, clusterNum, cluster2data, Pair_Index, &cn]() {
			//VectorXi seed;
			//MatrixXd C;
			//seed.resize(clusterNum);
			//for (int i = 0; i < clusterNum; i++)
			//{
			//	seed(i) = i;
			//}
			//jet(seed, true, C);

			//for (int i = 0; i < clusterNum; i++)
			//{
			//	VectorXi vertex_in_cluster;
			//	Extract_Vertex(cluster2data, Pair_Index, i, vertex_in_cluster);
			//	int num = vertex_in_cluster.rows();
			//	MatrixXd V_in_Cluster;
			//	V_in_Cluster.resize(num, 3);
			//	for (int j = 0; j < num; j++)
			//	{
			//		V_in_Cluster.row(j) = V.row(vertex_in_cluster(j));
			//	}
			//	m_viewer.data.add_points(V_in_Cluster, C.row(i));
			//}
			//m_viewer.data.points.resize(0, 0);
			m_viewer.data.lines.resize(0, 0);
			VectorXi vertex_in_cluster;
			Extract_Vertex(cluster2data, Pair_Index, cn, vertex_in_cluster);
			int num = vertex_in_cluster.rows();
			MatrixXd V_in_Cluster;
			V_in_Cluster.resize(num, 3);
			for (int j = 0; j < num; j++)
			{
				V_in_Cluster.row(j) = V.row(vertex_in_cluster(j));
			}
			m_viewer.data.set_points(V_in_Cluster, Eigen::RowVector3d(0, 1, 0));
			for (int j = 0; j < num/2; j++)
			{
				m_viewer.data.add_edges
				(
					V.row(vertex_in_cluster(j * 2)),
					V.row(vertex_in_cluster(j * 2 + 1)),
					Eigen::RowVector3d(1, 0, 0)
				);
			}
			/*m_viewer.data.add_edges
			(
				RowVector3d(0, 0, 0),
				RowVector3d(3, 0, 0),
				Eigen::RowVector3d(1, 0, 0)
			);
			m_viewer.data.add_edges
			(
				RowVector3d(0, 0, 0),
				RowVector3d(0, 3, 0),
				Eigen::RowVector3d(0, 1, 0)
			);
			m_viewer.data.add_edges
			(
				RowVector3d(0, 0, 0),
				RowVector3d(0, 0, 3),
				Eigen::RowVector3d(0, 0, 1)
			);*/
			cout << "clsuter " << cn << ": ";
			for (int j = 0; j < clusterCenter.cols(); j++)
			{
				cout << clusterCenter(cn, j) << ", ";
			}
			cout << endl;
			cn++;
			if (cn >= clusterNum)
			{
				cn = 0;
			}
		});

		// Generate menu
		m_viewer.screen->performLayout();

		return false;
	};

	m_viewer.launch();
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
	int percentage = 2;

	int pruned_point_size = V_Index_Prune.rows();
	int point_size = SV_col.rows();
	int p_dot_size = pruned_point_size / percentage;
	floor_index.resizeLike(V_Index_Prune);
	celling_index.resizeLike(V_Index_Prune);
	for (int i = 0; i < pruned_point_size; i++)
	{
		floor_index(i) = -1;
		celling_index(i) = -1;
	}

	//TODO: subset of P!
	for (int j = 0; j < p_dot_size; j++)
	{
		int i = j * percentage;
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

		for (int p = 0; p < point_size; p++)
		{
			double _sv_current = SV_Sorted(p);
			if (_sv_current >= _sv_floor)
			{
				floor_index(i) = p;
				break;
			}
		}
		for (int q = point_size - 1; q >= 0; q--)
		{
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

	//if (dimension == 7)
	//{
	//	VectorXi floor_index, celling_index;
	//	find_similar(V_Index_Prune, SV.col(0), SV_Sorted_7d, pair_threshold, floor_index, celling_index);
	//	for (int i = 0; i < pruned_point_size; i++)
	//	{
	//		for (int k = floor_index(i); k <= celling_index(i); k++)
	//		{
	//			if (V_Index_Dic(V_Index_Sorted_7d(k)) != -1 && V_Index_Sorted_7d(k) != V_Index_Prune(i))
	//			{
	//				_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted_7d(k);
	//			}
	//		}
	//	}
	//}
	//else if (dimension == 6)
	//{
	//	VectorXi floor_index_PV1, floor_index_PV2, celling_index_PV1, celling_index_PV2;
	//	find_similar(V_Index_Prune, SV.col(1), SV_Sorted_PV1, pair_threshold, floor_index_PV1, celling_index_PV1);
	//	find_similar(V_Index_Prune, SV.col(2), SV_Sorted_PV2, pair_threshold, floor_index_PV2, celling_index_PV2);
	//	for (int i = 0; i < pruned_point_size; i++)
	//	{	
	//		for (int m = floor_index_PV1(i); m <= celling_index_PV1(i); m++)
	//		{
	//			int V_Index_Sorted_PV1_ = V_Index_Sorted_PV1(m);
	//			for (int n = floor_index_PV2(i); n <= celling_index_PV2(i); n++)
	//			{
	//				int V_Index_Sorted_PV2_ = V_Index_Sorted_PV2(n);
	//				if (V_Index_Sorted_PV1_ == V_Index_Sorted_PV2_)
	//				{
	//					if (V_Index_Dic(V_Index_Sorted_PV1_) != -1 && V_Index_Sorted_PV1_ != V_Index_Prune(i))
	//					{
	//						//cout << SV(V_Index_Prune(i), 1) << "," << SV(V_Index_Prune(i), 2) << " and "<<SV(V_Index_Sorted_PV1_, 1) << "," << SV(V_Index_Sorted_PV1_, 2) << endl;
	//						_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted_PV1_;
	//					}
	//					break;
	//				}
	//			}
	//		}
	//	}
	//}
	if (dimension == 7)
	{
		VectorXi floor_index, celling_index;
		find_similar(V_Index_Prune, SV.col(0), SV_Sorted_7d, pair_threshold, floor_index, celling_index);
		for (int i = 0; i < pruned_point_size; i++)
		{
			if (floor_index(i) != -1 && celling_index(i) != -1)
			{
				for (int k = floor_index(i); k <= celling_index(i); k++)
				{
					if (V_Index_Dic(V_Index_Sorted_7d(k)) != -1 && V_Index_Sorted_7d(k) != V_Index_Prune(i))
					{
						_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted_7d(k);
					}
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
			if (floor_index_PV1(i) != -1 && celling_index_PV1(i) != -1)
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
	}
	//cout << "pair_size = " << pair_size << endl;
	Pair_Index = _Pair_Index.block(0, 0, pair_size, 2);
}

void Show_Pair_Points(const MatrixXd& V, const MatrixXi& Pair_Index, int num)
{
	if (num == -1)
	{
		for (int i = 0; i < Pair_Index.rows(); i++)
		{
			m_viewer.data.add_edges
			(
				V.row(Pair_Index(i, 0)),
				V.row(Pair_Index(i, 1)),
				Eigen::RowVector3d(1, 0, 0)
			);
		}
	}
	else
	{
		{
			m_viewer.data.lines.resize(0, 0);
			m_viewer.data.add_edges
			(
				V.row(Pair_Index(num, 0)),
				V.row(Pair_Index(num, 1)),
				Eigen::RowVector3d(1, 0, 0)
			);
		}
	}
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

void compute_rotation(const VectorXd& V_AXIS_start, const VectorXd& V_AXIS_end, const double angle, MatrixXd& Rot, MatrixXd& quaternionRot)
{
	//MatrixXd quaternionRot;
	quaternionRot.setIdentity(4, 4);
	VectorXd axisVect = V_AXIS_end - V_AXIS_start;
	axisVect.normalize();
	double cosA = cos(angle);
	double oneC = 1 - cosA;
	double sinA = sin(angle);
	double ux = axisVect(0);
	double uy = axisVect(1);
	double uz = axisVect(2);
	quaternionRot(0, 0) = ux*ux*oneC + cosA;
	quaternionRot(0, 1) = ux*uy*oneC - uz*sinA;
	quaternionRot(0, 2) = ux*uz*oneC + uy*sinA;
	quaternionRot(1, 0) = uy*ux*oneC + uz*sinA;
	quaternionRot(1, 1) = uy*uy*oneC + cosA;
	quaternionRot(1, 2) = uy*uz*oneC - ux*sinA;
	quaternionRot(2, 0) = uz*ux*oneC - uy*sinA;
	quaternionRot(2, 1) = uz*uy*oneC + ux*sinA;
	quaternionRot(2, 2) = uz*uz*oneC + cosA;

	
	//Transformation to origin
	MatrixXd _Trans2Origin;
	_Trans2Origin.setIdentity(4, 4);
	_Trans2Origin(0, 3) = -V_AXIS_start(0);
	_Trans2Origin(1, 3) = -V_AXIS_start(1);
	_Trans2Origin(2, 3) = -V_AXIS_start(2);

	MatrixXd _Trans2Axis;
	_Trans2Axis.setIdentity(4, 4);
	_Trans2Axis(0, 3) = V_AXIS_start(0);
	_Trans2Axis(1, 3) = V_AXIS_start(1);
	_Trans2Axis(2, 3) = V_AXIS_start(2);

	MatrixXd _Rot;
	_Rot.resize(4, 4);
	_Rot = _Trans2Axis * quaternionRot * _Trans2Origin;

	//Rot = _Rot.block(0, 0, 3, 3);
	Rot = _Rot;
	/*cout << "Rot = ";*/
	//printMatrix(Rot);
	//Vector4d zero;
	//zero << V_AXIS_start(0), V_AXIS_start(1), V_AXIS_start(2), 1;
	//zero = _Trans2Origin * zero;
	//cout << "test result 1:" << zero(0) << "," << zero(1) << "," << zero(2) << endl;
	//zero = quaternionRot * zero;
	//cout << "test result 1:" << zero(0) << "," << zero(1) << "," << zero(2) << endl;
	//zero = _Trans2Axis * zero;
	//cout << "test result 1:" << zero(0) << "," << zero(1) << "," << zero(2) << endl;
}

void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int dimension, MatrixXd& Transformation)
{
	int size = Pair_Index.rows();
	Transformation.resize(size, dimension);
	
	for (int i = 0; i < size; i++)
	{
		MatrixXd Frame_p1, Frame_p2;
		MatrixXd Rot, Rot1, quaternionRot_1, Rot2, quaternionRot_2;
		VectorXd Eular;
		VectorXd Trans;
		double Scale;
		Frame_p1.resize(3, 3);
		Frame_p2.resize(3, 3);
		Rot.resize(4, 4);
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

		//Frame_p1.col(0) = LN.row(pair_p1).transpose();
		//Frame_p1.col(1) = PD1.row(pair_p1).transpose();
		//Frame_p1.col(2) = PD2.row(pair_p1).transpose();

		//Frame_p2.col(0) = LN.row(pair_p2).transpose();
		//Frame_p2.col(1) = PD1.row(pair_p2).transpose();
		//Frame_p2.col(2) = PD2.row(pair_p2).transpose();

		//Rot = Frame_p2 * Frame_p1.inverse();
		//Trans = V.row(pair_p2).transpose() - Scale * Rot * (V.row(pair_p1).transpose());

		//same as paper
		VectorXd p1_p2_n;
		double cos_n1_n2, cos_p1_p2;
		double angle_n1_n2, angle_p1_p2;

		Vector3d N1 = LN.row(pair_p1), N2 = LN.row(pair_p2);
		p1_p2_n = N1.cross(N2);
		p1_p2_n.normalize();
		p1_p2_n = p1_p2_n + (Vector3d)V.row(pair_p1);
		cos_n1_n2 = (LN.row(pair_p1) * LN.row(pair_p2).transpose())(0,0) / (LN.row(pair_p1).norm()*LN.row(pair_p2).norm());
		angle_n1_n2 = acos(cos_n1_n2);
		compute_rotation(V.row(pair_p1), p1_p2_n, angle_n1_n2, Rot1, quaternionRot_1);


		VectorXd P_rot, P_N, P_D1, N_rot, P_D1_rot;
		N_rot.resize(4);
		P_rot.resize(4);
		P_D1.resize(4);
		P_N.resize(4);
		P_D1_rot.resize(4);
		P_rot << V(pair_p1, 0), V(pair_p1, 1), V(pair_p1, 2), 1;
		N_rot << LN(pair_p1, 0), LN(pair_p1, 1), LN(pair_p1, 2), 1;
		P_D1 << PD1(pair_p1, 0), PD1(pair_p1, 1), PD1(pair_p1, 2), 1;
		P_N = P_rot + N_rot;
		P_D1 = P_rot + P_D1;
		P_N(3) = 1;
		P_D1(3) =1;
		P_N = Rot1*P_N;
		P_D1 = Rot1 * P_D1;
		N_rot = P_N - P_rot;
		P_D1_rot = P_D1 - P_rot;

		MatrixXd  P_D1_rot_;
		P_D1_rot_.resize(1, 3);
		P_D1_rot_.row(0) << P_D1_rot(0), P_D1_rot(1), P_D1_rot(2);
		cos_p1_p2 = (P_D1_rot_.row(0) * (PD1.row(pair_p2)).transpose())(0,0) / (P_D1_rot_.norm() * (PD1.row(pair_p2)).norm());
		angle_p1_p2 = acos(cos_p1_p2);
		VectorXd normal = ((Vector3d)P_D1_rot_.row(0)).cross((Vector3d)PD1.row(pair_p2));
		//cout << "normal: ";
		//printVector(normal);
		//cout << "N_ROT: ";
		//printVector(N_rot);
		if (normal.dot(N_rot.head(3)) < 0)
		{
			//cout << "different" << endl;
			angle_p1_p2 = -angle_p1_p2;
		}
		compute_rotation(V.row(pair_p1), V.row(pair_p1) + N_rot.head(3).transpose(), angle_p1_p2, Rot2, quaternionRot_2);
		//Rot = Rot2 * Rot1;
		Rot = quaternionRot_2.block(0, 0, 3, 3) * quaternionRot_1.block(0, 0, 3, 3);
		Trans = V.row(pair_p2).transpose() - Scale * (V.row(pair_p1).transpose());

		

		//VectorXd test_point;
		//test_point.resize(4);
		//test_point << V(pair_p1, 0), V(pair_p1, 1), V(pair_p1, 2), 1;
		//cout << "before point:" << endl;
		//printVector(test_point);
		//test_point = Rot*test_point;
		//cout << "after point:" << endl;
		//printVector(test_point);

		//cout << "Rot BEFORE= " << endl
		//	<< Rot(0, 0) << "," << Rot(0, 1) << "," << Rot(0, 2) << "," << endl
		//	<< Rot(1, 0) << "," << Rot(1, 1) << "," << Rot(1, 2) << "," << endl
		//	<< Rot(2, 0) << "," << Rot(2, 1) << "," << Rot(2, 2) << "," << endl;

		RotateMatrix2EularAngle(Rot, Eular);
		
		/*MatrixXd Rot_TEST;
		
		RotateEularAngle2Matrix(Eular, Rot_TEST);
		cout << "Rot after= " << endl
			<< Rot_TEST(0, 0) << "," << Rot_TEST(0, 1) << "," << Rot_TEST(0, 2) << "," << endl
			<< Rot_TEST(1, 0) << "," << Rot_TEST(1, 1) << "," << Rot_TEST(1, 2) << "," << endl
			<< Rot_TEST(2, 0) << "," << Rot_TEST(2, 1) << "," << Rot_TEST(2, 2) << "," << endl;*/

	/*	VectorXd P_D1_, P_D1_rot_ROT;
		P_D1_.resize(4);
		P_D1_rot_ROT.resize(4);
		P_D1_ = Rot2 * P_D1;
		P_D1_rot_ROT = P_D1_ - P_rot;*/

	/*	cout << "Rot AIM= " << endl
			<< LN.row(pair_p2)(0) << "," << LN.row(pair_p2)(1) << "," << LN.row(pair_p2)(2) << "," << endl
			<< PD1.row(pair_p2)(0) << "," << PD1.row(pair_p2)(1) << "," << PD1.row(pair_p2)(2) << "," << endl
			<< V(pair_p2, 0) << "," << V(pair_p2, 1) << "," << V(pair_p2, 2) << "," << endl;*/

		/*cout << "Rot BEFORE= " << endl
			<< LN.row(pair_p1)(0) << "," << LN.row(pair_p1)(1) << "," << LN.row(pair_p1)(2) << "," << endl
			<< PD1.row(pair_p1)(0) << "," << PD1.row(pair_p1)(1) << "," << PD1.row(pair_p1)(2) << "," << endl
			<< V(pair_p1, 0) << "," << V(pair_p1, 1) << "," << V(pair_p1, 2) << "," << endl;*/

		//cout << "Rot AFTER= " << endl
		//	<< N_rot(0) << "," << N_rot(1) << "," << N_rot(2) << "," << N_rot(3) << endl
		//	<< P_D1_rot_ROT(0) << "," << P_D1_rot_ROT(1) << "," << P_D1_rot_ROT(2) << "," << P_D1_rot_ROT(3) << endl;
			//<< N_rot(0) << "," << N_rot(1) << "," << N_rot(2) << "," << N_rot(3) << endl;

		
		//Trans = V.row(pair_p2).transpose() - Scale * Rot.block(0,0,3,3) * (V.row(pair_p1).transpose());

		//cout <<"V1="<< V(pair_p1, 0) << "," << V(pair_p1, 1) << "," << V(pair_p1, 2) << "," << endl;
		//cout <<"V2="<< V(pair_p2, 0) << "," << V(pair_p2, 1) << "," << V(pair_p2, 2) << "," << endl;
		
		//cout << "TRANS=" << Trans(0) << "," << Trans(1) << "," << Trans(2) << "," << endl;

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

void RotateMatrix2EularAngle(const MatrixXd& Rot, VectorXd& Eular)
{
	Eular.resize(3);
	double 
	r1 = Rot(0, 0),  r2 = Rot(0, 1),   r3 = Rot(0, 2),
	r4 = Rot(1, 0),  r5 = Rot(1, 1),   r6 = Rot(1, 2),
	r7 = Rot(2, 0),  r8 = Rot(2, 1),   r9 = Rot(2, 2);
	double Rx = atan2(r4 , r1);
	double sx = sin(Rx), cz = cos(Rx);
	double Ry = atan2((-r7) , (r1*cz + r4*sx));
	double Rz = atan2((r3*sx - r6*cz) , (r5*cz - r2*sx));
	//cout << Rx << "," << Ry << "," << Rz << endl;
	Eular << Rx, Ry, Rz;
}

void RotateEularAngle2Matrix(const VectorXd& Eular, MatrixXd& Rot)
{
	Rot.resize(3, 3);
	double E[3];

	int i = 0;
	for (i = 0; i<3; i++)
	{
		E[i] = Eular(i);
	}
	/*N*/
	Rot(0, 0) = 1.0 * (cos(E[0])*cos(E[1]));
	Rot(1, 0) = 1.0 * (sin(E[0])*cos(E[1]));
	Rot(2, 0) = -1.0 * (sin(E[1]));
	/*O*/
	Rot(0, 1) = 1.0 * (cos(E[0])*sin(E[1])*sin(E[2]) - sin(E[0])*cos(E[2]));
	Rot(1, 1) = 1.0 * (sin(E[0])*sin(E[1])*sin(E[2]) + cos(E[0])*cos(E[2]));
	Rot(2, 1) = 1.0 * (cos(E[1])*sin(E[2]));
	/*A*/
	Rot(0, 2) = 1.0 * (cos(E[0])*sin(E[1])*cos(E[2]) + sin(E[0])*sin(E[2]));
	Rot(1, 2) = 1.0 * (sin(E[0])*sin(E[1])*cos(E[2]) - cos(E[0])*sin(E[2]));
	Rot(2, 2) = 1.0 * (cos(E[1])*cos(E[2]));

}

int MeanShift_Cluster(const MatrixXd& V, const MatrixXd& dataPts, MatrixXd& clusterCenter, MatrixXd& data2cluster, MatrixXd& cluster2data)
{
	// Matlab instance
	Engine* engine;

	// Launch MATLAB
	matlab::mlinit(&engine);

	// Send V matrix to matlab for bounding box size
	MatrixXd V_ = V.transpose();
	matlab::mlsetmatrix(&engine, "V", V_);
	matlab::mleval(&engine, "[boxSize] = BoundingBoxSize(V);");
	cout << "box size = " << matlab::mlgetscalar(&engine, "boxSize") << endl;

	// Send Transformation matrix to matlab
	MatrixXd dataPts_ = dataPts.transpose();
	matlab::mlsetmatrix(&engine, "dataPts", dataPts_);
	
	//Cluster using matlab
	matlab::mleval(&engine, "[clustCent,data2cluster,cluster2data,numClust] = MeanShiftCluster(dataPts,1,boxSize);");
	
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
	
	// Close MATLAB
	matlab::mlclose(&engine);

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

void printVector(const VectorXi& vec)
{
	int size = vec.rows();
	for (int i = 0; i < size; i++)
	{
		cout << vec(i) << ",";
	}
	cout << endl;
}

void printVector(const VectorXd& vec)
{
	int size = vec.rows();
	for (int i = 0; i < size; i++)
	{
		cout << vec(i) << ",";
	}
	cout << endl;
}

void printMatrix(const MatrixXd& mat)
{
	int r = mat.rows(), c = mat.cols();
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c - 1; j++)
		{
			cout << mat(i, j) << ",";
		}
		cout << mat(i, c - 1) << endl;
	}
}

void printMatrix(const MatrixXi& mat)
{
	int r = mat.rows(), c = mat.cols();
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c - 1; j++)
		{
			cout << mat(i, j) << ",";
		}
		cout << mat(i, c - 1) << endl;
	}
}