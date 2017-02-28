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
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>
#include "symDetec_shared_path.h"
using namespace igl;
using namespace Eigen;
using namespace std;

void RotateMatrix2EularAngle(const MatrixXd Rot, VectorXd& Eular);
void prune_umbilic_points(const VectorXd& SV, const VectorXi& V_Index_Sample, VectorXi& V_Index_Prune, VectorXi& V_Index_Dic, double prune_threshold);
void compute_signature(const VectorXd& PV1, const VectorXd& PV2, VectorXd& SV);
void samping(int size, double Sample_Rate, VectorXi& V_Index_Sample);
void search_pair(const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted, const VectorXd& SV, const VectorXd& SV_Sorted, double pair_threshold, MatrixXi& Pair_Index);
void compute_local_frame(const MatrixXd& PD1, const MatrixXd& PD2, MatrixXd& LN);
void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int dimension, MatrixXd& Transformation);
int MeanShift_Cluster(const MatrixXd& dataPts, double bandwidth, MatrixXd& clusterCenter, MatrixXi& data2cluster, MatrixXi& cluster2dataCell);

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
	VectorXd SV, SV_Sorted;

	VectorXi V_Index_Sample, V_Index_Dic, V_Index_Prune, V_Index_Sorted;
	MatrixXi Pair_Index;

	MatrixXd Transformation;

	int clusterNum;
	MatrixXd clusterCenter;
	MatrixXi data2cluster, cluster2dataCell;

	// Load a mesh in OFF format
	readOBJ(SYMDETEC_SHARED_PATH "/bunny.OBJ", V, F);
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
	samping(V.rows(), 0.2, V_Index_Sample);
	cout << "SAMPLE FINISH!" << endl;
	cout << "After Sampling, Vertex: " << V_Index_Sample.rows() << endl << endl;

	//Prune umbilic points
	prune_umbilic_points(SV, V_Index_Sample, V_Index_Prune, V_Index_Dic, prune_threshold);
	cout << "PRUNE FINISH!" << endl;
	cout << "After Pruning, Vertex: " << V_Index_Prune.rows() << endl << endl;

	sortrows(SV, true, SV_Sorted, V_Index_Sorted);
	//for (int i = 0; i < SV.rows(); i++)
	//{
	//	cout << SV_Sorted(i) <<" , " << SV(V_Index_Sorted(i)) << endl;
	//}

	//Pairing
	search_pair(V_Index_Dic, V_Index_Prune, V_Index_Sorted, SV, SV_Sorted, pair_threshold, Pair_Index);
	cout << "PAIR FINISH!" << endl;
	cout << "After Pairing, Pair: " << Pair_Index.rows() << endl << endl;

	//Compute Transformation 
	compute_transformation(Pair_Index, V, PV1, PV2, PD1, PD2, LN, dimension, Transformation);
	cout << "COMPUTE TRANSFORMATION FINISH!" << endl << endl;

	//Mean-Shift
	clusterNum = MeanShift_Cluster(Transformation, bandwidth, clusterCenter, data2cluster, cluster2dataCell);
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

	//for (int i = 0; i < data2cluster.rows(); i++)
	//{
	//	cout << "point " << i << " ---> " << data2cluster(i,0) << endl;
	//}
	cout << endl;


	//// Init the viewer
	//viewer::Viewer viewer;
	//// Plot the mesh
	//viewer.data.set_mesh(V, F);
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
	//viewer.launch();
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

void prune_umbilic_points(const VectorXd& SV, const VectorXi& V_Index_Sample, VectorXi& V_Index_Prune, VectorXi& V_Index_Dic, double prune_threshold)
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
		if (fabs(SV(V_Index_Sample(i))) < prune_threshold)
		{
			V_Index_Dic(V_Index_Sample(i)) = size_P;
			_V_Index_Prune(size_P++) = V_Index_Sample(i);
			//cout << "add " << i << endl;
		}
	}

	V_Index_Prune = _V_Index_Prune.head(size_P);
}

void compute_signature(const VectorXd& PV1, const VectorXd& PV2, VectorXd& SV)
{
	int size = PV1.rows();
	SV.resizeLike(PV1);

	for (int i = 0; i < size; i++)
	{
		if(PV2(i) == 0.0)
		{
			SV(i) = (PV1(i) > 0.0) ? (-INF) : (INF);
		}
		else
		{
			SV(i) = PV1(i) / PV2(i);
		}

		//cout << "SV" << i << " = " << SV(i) << endl;
	}

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

void search_pair(const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted, const VectorXd& SV, const VectorXd& SV_Sorted, double pair_threshold, MatrixXi& Pair_Index)
{
	int point_size = SV.rows();
	int pruned_point_size = V_Index_Prune.rows();
	int max_size = pruned_point_size * (pruned_point_size - 1) / 2;
	int pair_size = 0;

	MatrixXi _Pair_Index(max_size, 2);

	for (int i = 0; i < pruned_point_size; i++)
	{
		double _sv = SV(V_Index_Prune(i));
		//cout << "sv = " << _sv << endl;
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

		//cout << "floor = " << _sv_floor << " celling = " << _sv_celling << endl;
		int floor_index = -1, celling_index = -1;
		for (int p = 0; p < point_size; p++) {
			double _sv_current = SV_Sorted(p);
			if (_sv_current >= _sv_floor)
			{
				//cout << "find floor = " << _sv_current << endl;
				floor_index = p;
				break;
			}
		}
		for (int q = point_size - 1; q >= 0; q--) {
			double _sv_current = SV_Sorted(q);
			if (_sv_current <= _sv_celling)
			{
				//cout << "find celling = " << _sv_current << endl;
				celling_index = q;
				break;
			}
		}
		//cout << floor_index << " , " << celling_index << endl;
		for (int k = floor_index; k < celling_index; k++)
		{
			//cout << "SV" << V_Index_Prune(i) << " = " << SV(V_Index_Prune(i)) << endl;
			//cout << "SV" << V_Index_Sorted(k) << " = " << SV(V_Index_Sorted(k)) << endl;
			//cout << "DIC" << V_Index_Sorted(k) << " = " << V_Index_Dic(V_Index_Sorted(k)) << endl;
			if (V_Index_Dic(V_Index_Sorted(k)) != -1)
			{
				//cout << "add pair" << V_Index_Prune(i) << " , " << V_Index_Sorted(k) << endl;
				_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted(k);
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

int MeanShift_Cluster(const MatrixXd& dataPts, double bandwidth, MatrixXd& clusterCenter, MatrixXi& data2cluster, MatrixXi& cluster2dataCell)
//void MeanShift_Cluster()
{
//	int numDim = dataPts.cols();
//	int numPts = dataPts.rows();
//	int numCluster = 0;
//	int numInitPts = numPts;
//	double bandSq = bandwidth * bandwidth;
//	double stopThreshold = 1e-3 * bandwidth;
//	VectorXi visitedFlag, clusterVotes, initPtIndex;
//	visitedFlag.setZero(numPts);
//	clusterVotes.setZero(numPts);
//	initPtIndex.resize(numPts);
//	for (int i = 0; i < numPts; i++)
//	{
//		initPtIndex(i) = i;
//	}
//
//	while (numInitPts)
//	{
//		int randomSeedIndex = ceil((numInitPts - 1e-6) * rand());
//		int startIndex = initPtIndex(randomSeedIndex);
//		VectorXd myMean = dataPts.row(startIndex);
//		VectorXi myMembers, thisClusterVotes;
//		myMembers.resize(numPts);
//		thisClusterVotes.setZero(numPts);
//
//		while (1)
//		{
//			VectorXd sqDistToAll; //dist squared from mean to all points still active
//			VectorXd inInds = find(sqDistToAll < bandSq); //points within bandWidth
//		
//		}
//
//	}
	// Matlab instance
	Engine* engine;
	//igl::MatlabWorkspace mw;
	//mw.save(V, "V");
	//mw.save_index(F, "F");
	//mw.save(L, "L");
	//mw.write("fertility.mat");

	// Launch MATLAB
	matlab::mlinit(&engine);

	// Send Transformation matrix to matlab
	MatrixXd dataPts_ = dataPts.transpose();
	matlab::mlsetmatrix(&engine, "dataPts", dataPts_);

	// Send bandwidth to matlab
	matlab::mlsetscalar(&engine, "bandwidth", bandwidth);

	//Cluster using matlab
	igl::matlab::mleval(&engine, "[data2cluster,cluster2dataCell,clustCent,numClust] = MeanShiftCluster(dataPts,bandwidth)");
	
	// Get clustering
	matlab::mlgetmatrix(&engine, "clustCent", clusterCenter);
	matlab::mlgetmatrix(&engine, "data2cluster", data2cluster);
	matlab::mlgetmatrix(&engine, "cluster2dataCell", cluster2dataCell);
	
	int numClust = -1;
	numClust = (int)matlab::mlgetscalar(&engine, "numClust");
	//cout << "matlab" << matlab::mleval(&engine, "size(clustCent)") << endl;
	//cout << "matlab" << matlab::mleval(&engine, "size(data2cluster)") << endl;
	//cout << "matlab" << matlab::mleval(&engine, "size(cluster2dataCell)") << endl;
	//cout << data2cluster.rows() << " , " << data2cluster.cols() << endl;

	//cout << cluster2dataCell.rows() << " , " << cluster2dataCell.cols() << endl;
	clusterCenter.transposeInPlace();
	data2cluster.transposeInPlace();
	return numClust;
}