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
#include <igl/doublearea.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <iostream>


using namespace igl;
using namespace Eigen;
using namespace std;

#define TRANSFORM_7D		     0
#define TRANSFORM_6D	         1
#define ROTATION_3D	             2
#define SYMMETRY_PLANE_4D        3
#define ROTATION_AND_SYMMETRY_7D 4
#define SYMMETRY_VECTOR_3D       5
#define ROTATION_AND_SVECTOR_6D  6

//// Init the viewer
viewer::Viewer m_viewer;

// Matlab instance
Engine* engine;

void RotateMatrix2EularAngle(const MatrixXd& Rot, VectorXd& Eular);
void RotateEularAngle2Matrix(const VectorXd& Eular, MatrixXd& Rot);
void prune_umbilic_points(const MatrixXd& SV, const VectorXi& V_Index_Sample, VectorXi& V_Index_Prune, VectorXi& V_Index_Dic, double prune_threshold);
double compute_boxSize(const MatrixXd& V);
void compute_signature(const VectorXd& PV1, const VectorXd& PV2, MatrixXd& SV);
void random_sample(int size, double Sample_Rate, VectorXi& V_Index_Sample);
void sort_vertex_in_face_by_row(MatrixXi& F);
void sort_vertex_in_face_by_col(int start, int end, int col, MatrixXi& F);
void sort_vertex_in_face(int vertexSize, MatrixXi& F);
void sample_based_on_area(double ratio, const MatrixXd& V, MatrixXi& F, VectorXi& V_Index_Sample);
void sampleFromFaceGroup(int samplePerUnit, const MatrixXi& F, int SampleFaceStart, int SampleFaceEnd, int BigFaceNum, const VectorXi& BigFaceIndex, int& sampleNum, VectorXi& V_Sample);
void sort_by_signature(const MatrixXd& SV, VectorXd& SV_Sorted_7d, VectorXd& SV_Sorted_PV1, VectorXd& SV_Sorted_PV2, VectorXi& V_Index_Sorted_7d, VectorXi& V_Index_Sorted_PV1, VectorXi& V_Index_Sorted_PV2);
void find_similar(const VectorXi& V_Index_Prune, const VectorXd& SV_col, const VectorXd& SV_Sorted, double pair_threshold, VectorXi& floor_index, VectorXi& celling_index);
void search_pair(int mode, double boxSize, const MatrixXd& V, const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted_7d, const VectorXi& V_Index_Sorted_PV1, const VectorXi& V_Index_Sorted_PV2, const MatrixXd& SV, const VectorXd& SV_Sorted_7d, const VectorXd& SV_Sorted_PV1, const VectorXd& SV_Sorted_PV2, double pair_threshold, MatrixXi& Pair_Index);
void compute_curvature(const MatrixXd& V, const MatrixXi& F, double radius, MatrixXd& VPD1, MatrixXd&PD2, VectorXd& PV1, VectorXd& PV2);
void compute_local_frame(const MatrixXd& PD1, const MatrixXd& PD2, MatrixXd& LN);
void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int mode, MatrixXd& Transformation);
void compute_rotation(const VectorXd& V_AXIS_start, const VectorXd& V_AXIS_end, const double angle, MatrixXd& Rot, MatrixXd& quaternionRot);
int  MeanShift_Cluster(int mode,const MatrixXd& V, const MatrixXd& Transformation, MatrixXd& clusterCenter, MatrixXd& data2cluster, MatrixXd& cluster2data);
void patching(const MatrixXd& V, const MatrixXi& F, const MatrixXi& Pair_Index, const MatrixXd& clusterCenter, const MatrixXd& cluster2data, MatrixXd& Patch1, MatrixXd& Patch2, VectorXi& size, VectorXi& rank);
void Extract_Vertex(const MatrixXd& cluster2data, const MatrixXi& Pair_Index, int which_cluster, VectorXi& V_in_Cluster);
void Show_Pair_Points(const MatrixXd& V, const MatrixXi& Pair_Index, int num);
void showCluster(const MatrixXd& Transformation, const MatrixXd& clusterCenter, const MatrixXd& data2cluster);
void printVector(const VectorXi& vec);
void printVector(const VectorXd& vec);
void printVector(const Vector3d& vec);
void printMatrix(const MatrixXd& mat);
void printMatrix(const MatrixXi& mat);


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

	//cout << "size = " << size << endl;
	for (int i = 0; i < size_Sample; i++)
	{
		//cout << i << "i = " << SV(i) << endl;
		if (fabs(SV(V_Index_Sample(i),0)) < prune_threshold)
		{
			/*cout << "1 " <<V_Index_Sample(i) << endl;
			cout << "2 " <<V_Index_Dic(V_Index_Sample(i)) << endl;*/
			V_Index_Dic(V_Index_Sample(i)) = size_P;
			_V_Index_Prune(size_P++) = V_Index_Sample(i);
			//cout << "add " << i << endl;
		}
	}

	V_Index_Prune = _V_Index_Prune.head(size_P);
}

double compute_boxSize(const MatrixXd& V)
{
	// Send V matrix to matlab for bounding box size
	MatrixXd V_ = V.transpose();
	matlab::mlsetmatrix(&engine, "V", V_);
	matlab::mleval(&engine, "[boxSize] = BoundingBoxSize(V);");
	double boxSize = matlab::mlgetscalar(&engine, "boxSize");
	return boxSize;
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

void random_sample(int size, double Sample_Rate, VectorXi& V_Index_Sample)
{
	//just random sample part of points 
	int sample_size = size * Sample_Rate;
	V_Index_Sample.resize(sample_size);

	vector<int> s_stl;
	for (int i = 0; i<size; ++i) s_stl.push_back(i);
	random_shuffle(s_stl.begin(), s_stl.end());

	int* temp = new int[sample_size];

	for (int j = 0; j < sample_size; j++)
	{
		temp[j] = s_stl.at(j);
	}

	sort(temp, temp + sample_size);

	for (int j = 0; j < sample_size; j++)
	{
		V_Index_Sample(j) = temp[j];
	}
}
void sort_vertex_in_face_by_row(MatrixXi& F)
{
	int size = F.rows();
	for (int i = 0; i < size; i++)
	{
		vector<int> row, row_sorted;
		row.push_back(F(i, 0));
		row.push_back(F(i, 1));
		row.push_back(F(i, 2));
		row_sorted = row;
		sort(row_sorted.begin(), row_sorted.end());
		//cout << "before" << F.row(i) << endl;
		if (row_sorted[0] == row[0])
		{
			F.row(i) << row[0], row[1], row[2];
		}
		else if (row_sorted[0] == row[1])
		{
			F.row(i) << row[1], row[2], row[0];
		}
		else
		{
			F.row(i) << row[2], row[0], row[1];
		}
		//cout << "after" << F.row(i) << endl;
	}
		
}

void sort_vertex_in_face_by_col(int start, int end, int col, MatrixXi& F)
{

	MatrixXi F_ = F;
	VectorXi sortCol, sortedCol, sortedIndex;
	sortCol = F_.block(start, col, end - start + 1, 1);
	sortrows(sortCol, true, sortedCol, sortedIndex);
	for (int i = start; i <= end; i++)
	{
		F.row(i) = F_.row(start + sortedIndex(i - start));
	}
}

void sort_vertex_in_face(int vertexSize, MatrixXi& F)
{
	int size = F.rows();
	
	sort_vertex_in_face_by_row(F);

	sort_vertex_in_face_by_col(0, size - 1, 0, F);


	//int start, end, currentVertex;
	//start = 0;
	//currentVertex = F(0, 0);
	//for (int i = 1; i < F.rows(); i++)
	//{
	//	if (F(i, 0) != currentVertex)
	//	{
	//		end = i - 1;
	//		if (start < end)
	//		{
	//			sort_vertex_in_face_by_col(start, end, 1, F);
	//		}
	//		currentVertex = F(i, 0);
	//		start = i;
	//	}
	//}
	//printMatrix(F);
	

	//start = 0;
	//currentVertex = F(0, 1);
	//for (int i = 1; i < F.rows(); i++)
	//{
	//	if (F(i, 1) != currentVertex)
	//	{
	//		end = i - 1;
	//		if (start < end)
	//		{
	//			sort_vertex_in_face_by_col(start, end, 2, F);
	//		}
	//		currentVertex = F(i, 1);
	//		start = i;
	//	}
	//}
	//printMatrix(F);
	

}

void sample_based_on_area(double sampleRatio, const MatrixXd& V, MatrixXi& F, VectorXi& V_Index_Sample)
{
	VectorXi V_Sample;
	VectorXd DbArea;
	V_Sample.resize(V.rows());
	V_Sample.setZero();
	int sampleNum = 0;

	sort_vertex_in_face(V.rows(), F);

	doublearea(V, F, DbArea);

	VectorXd DbGroup_Area;
	DbGroup_Area.resize(V.rows());
	int currentVertex;
	currentVertex = F(0, 0);
	DbGroup_Area(0) = 0.0;
	double min = DbArea(0), max = 0.0;
	for (int i = 0; i < F.rows(); i++)
	{
		if (F(i, 0) == currentVertex)
		{
			DbGroup_Area(currentVertex) += DbArea(i);
		}
		else
		{
			if (DbGroup_Area(currentVertex) != 0 && min > DbGroup_Area(currentVertex))
			{
				min = DbGroup_Area(currentVertex);
			}
			if (max < DbGroup_Area(currentVertex))
			{
				max = DbGroup_Area(currentVertex);
			}
			currentVertex = F(i, 0);
			DbGroup_Area(currentVertex) = 0.0;
		}
	}

	//cout << "min = " << min << endl;
	//cout << "max = " << max << endl;

	double totalArea = 0.0;
	for (int i = 0; i < F.rows(); i++)
	{
		totalArea += DbArea(i);
	}
	//cout << "total = " << totalArea << endl;
	//cout << F.rows() << endl;
	//cout << totalArea / F.rows() << endl;
	//int temp = F.rows();
	//int width = 0;
	//while (temp / 10 > 1) {
	//	width++;
	//	temp /= 10;
	//}
	
	//int sampleUnits = pow(10,width - 1);
	//double unitArea = totalArea / sampleUnits;
	double unitArea = 10 * max;
	int sampleUnits = ceil(totalArea / unitArea);
	int samplePerUnit = ceil((double)V.rows() * sampleRatio / (double)sampleUnits);
	//cout << sampleUnits << endl;
	//cout << samplePerUnit << endl;
	//cout << "unitArea = " << unitArea << endl;

	double currentGroupArea;
	VectorXi BigFaceIndex;
	int BigFaceNum;
	BigFaceIndex.resize(V.rows());
	int currentGroupStart, currentGroupEnd;
	currentGroupStart = currentGroupEnd = -1;
	for (int i = 0; i < F.rows(); i++)
	{
		double DbAreai = DbArea(i);
		if (currentGroupStart < 0)
		{
			currentGroupArea = 0;
			currentGroupStart = currentGroupEnd = i;
			BigFaceNum = 0;
			BigFaceIndex.setOnes();
			BigFaceIndex = BigFaceIndex * -1;
		}

		if (DbAreai >= unitArea)
		{
			BigFaceIndex(BigFaceNum++) = i;
			//TODO extract from dbAI
			sampleFromFaceGroup(samplePerUnit, F, currentGroupStart, currentGroupEnd, BigFaceNum, BigFaceIndex, sampleNum, V_Sample);

		}
		else
		{
			currentGroupArea += DbArea(i);
			currentGroupEnd++;
			if (currentGroupArea >= unitArea)
			{
				//cout << "start = " << currentGroupStart << ", end = " << currentGroupEnd << endl;
				//TODO extract from start to end
				sampleFromFaceGroup(samplePerUnit, F, currentGroupStart, currentGroupEnd, BigFaceNum, BigFaceIndex, sampleNum, V_Sample);
				currentGroupStart = currentGroupEnd = -1;
			}
		}

	}

	//printVector(V_Sample);

	V_Index_Sample.resize(sampleNum);
	int current = 0;
	for (int i = 0; i < V.rows(); i++)
	{
		if (V_Sample(i) == 1)
		{
			V_Index_Sample(current++) = i;
		}
	}
}

void sampleFromFaceGroup(int samplePerUnit, const MatrixXi& F, int SampleFaceStart, int SampleFaceEnd, int BigFaceNum, const VectorXi& BigFaceIndex, int& sampleNum, VectorXi& V_Sample)
{
	int thisSampleNum;
	if (SampleFaceStart == SampleFaceEnd)
	{
		/*cout << "here" << endl;*/
		thisSampleNum = (samplePerUnit > 3) ? 3 : samplePerUnit;
		for (int i = 0; i < thisSampleNum; i++)
		{
			int thisPointIndex = F(SampleFaceStart, i);
			if (V_Sample(thisPointIndex) == 0) {
				V_Sample(thisPointIndex) = 1;
				sampleNum++;
			}
		}
	}
	else
	{
		vector<int> temp;
		//int* no_repeat = new int[3 * (SampleFaceEnd - SampleFaceStart)];
		int currentBigFace = 0;
		for (int i = SampleFaceStart; i <= SampleFaceEnd; i++)
		{
			if (BigFaceIndex(currentBigFace) != i)
			{
				for (int j = 0; j < 3; j++)
				{
					int thisPointIndex = F(i, j);
					if (V_Sample(thisPointIndex) == 0) {
						//cout << "maybe repeat " << thisPointIndex << endl;
						temp.push_back(thisPointIndex);
					}
				}
			}
			else
			{
				currentBigFace++;
			}
			
		}

		//delete repetation
		vector<int> no_repeat;
		sort(temp.begin(), temp.end());
		if (temp.size() == 0)
		{
			return;
		}
		int last = temp[0];
		for (int i = 1; i < temp.size(); i++)
		{
			if (last != temp[i])
			{
				no_repeat.push_back(last);
				last = temp[i];
			}
		}
		if (no_repeat.size() == 0 || no_repeat[no_repeat.size() - 1] != last)
		{
			no_repeat.push_back(last);
		}


		//random shuffle
		random_shuffle(no_repeat.begin(), no_repeat.end());
		
		//select sample and mark which to be sampled
		if (samplePerUnit > no_repeat.size())
		{
			thisSampleNum = no_repeat.size();
		}
		else
		{
			thisSampleNum = samplePerUnit;
		}
		//cout << "samplePerUnit = " << samplePerUnit << endl;
		//cout << "thisGroupSize = " << no_repeat.size() << endl;
		//cout << "thisSampleNum = " << thisSampleNum << endl;
		for (int i = 0; i < thisSampleNum; i++)
		{
			if (V_Sample(no_repeat[i]) == 0)
			{
				V_Sample(no_repeat[i]) = 1;
				sampleNum++;
			}
		}

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

void search_pair(int mode, double boxSize,const MatrixXd& V, const VectorXi& V_Index_Dic, const VectorXi& V_Index_Prune, const VectorXi& V_Index_Sorted_7d, const VectorXi& V_Index_Sorted_PV1, const VectorXi& V_Index_Sorted_PV2, const MatrixXd& SV, const VectorXd& SV_Sorted_7d, const VectorXd& SV_Sorted_PV1, const VectorXd& SV_Sorted_PV2, double pair_threshold, MatrixXi& Pair_Index)
{
	int point_size = SV.rows();
	int pruned_point_size = V_Index_Prune.rows();
	int max_size = pruned_point_size * (pruned_point_size - 1) / 2;
	int pair_size = 0;

	double dist_threshold = boxSize * 0.1;

	MatrixXi _Pair_Index(max_size, 2);

	if (mode == TRANSFORM_7D)
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
	else
	{
		VectorXi floor_index_PV1, floor_index_PV2, celling_index_PV1, celling_index_PV2;
		find_similar(V_Index_Prune, SV.col(1), SV_Sorted_PV1, pair_threshold, floor_index_PV1, celling_index_PV1);
		//find_similar(V_Index_Prune, SV.col(2), SV_Sorted_PV2, pair_threshold, floor_index_PV2, celling_index_PV2);
		for (int i = 0; i < pruned_point_size; i++)
		{
			int my_Index = V_Index_Prune(i);
			double this_point_pv1 = SV(my_Index, 1);
			double this_point_pv2 = SV(my_Index, 2);
			double pv1_min = this_point_pv1 * (1 - pair_threshold);
			double pv1_max = this_point_pv1 * (1 + pair_threshold);
			double pv2_min = this_point_pv2 * (1 - pair_threshold);
			double pv2_max = this_point_pv2 * (1 + pair_threshold);
			if (floor_index_PV1(i) != -1 && celling_index_PV1(i) != -1)
			{
				for (int m = floor_index_PV1(i); m <= celling_index_PV1(i); m++)
				{
					int V_Index_Sorted_PV1_ = V_Index_Sorted_PV1(m);
					double pv1 = SV(V_Index_Sorted_PV1_, 1);
					double pv2 = SV(V_Index_Sorted_PV1_, 2);
					if (pv2 >= pv2_min && pv2 <= pv2_max)
					{
						if (V_Index_Dic(V_Index_Sorted_PV1_) != -1 && V_Index_Sorted_PV1_ > my_Index)
						{
							Vector3d dist = V.row(my_Index) - V.row(V_Index_Sorted_PV1_);
							if (dist.norm() >= dist_threshold)
							{
								//cout << SV(V_Index_Prune(i), 1) << "," << SV(V_Index_Prune(i), 2) << " and "<<SV(V_Index_Sorted_PV1_, 1) << "," << SV(V_Index_Sorted_PV1_, 2) << endl;
								_Pair_Index.row(pair_size++) << my_Index, V_Index_Sorted_PV1_;
							}
							
						}
					}
					//cout << "pv1: " << pv1_min << " , " << pv1_max << " , " << pv1 << endl;
					//cout << "pv2: " << pv2_min << " , " << pv2_max << " , " << pv2 << endl;
					//for (int n = floor_index_PV2(i); n <= celling_index_PV2(i); n++)
					//{
					//	int V_Index_Sorted_PV2_ = V_Index_Sorted_PV2(n);
					//	if (V_Index_Sorted_PV1_ == V_Index_Sorted_PV2_)
					//	{
					//		if (V_Index_Dic(V_Index_Sorted_PV1_) != -1 && V_Index_Sorted_PV1_ != V_Index_Prune(i))
					//		{
					//			//cout << SV(V_Index_Prune(i), 1) << "," << SV(V_Index_Prune(i), 2) << " and "<<SV(V_Index_Sorted_PV1_, 1) << "," << SV(V_Index_Sorted_PV1_, 2) << endl;
					//			_Pair_Index.row(pair_size++) << V_Index_Prune(i), V_Index_Sorted_PV1_;
					//		}
					//		break;
					//	}
					//}
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
	
void compute_curvature(const MatrixXd& V, const MatrixXi& F, double radius, MatrixXd& PD1, MatrixXd&PD2, VectorXd& PV1, VectorXd& PV2)
{
	MatrixXd V_ = V.transpose();
	MatrixXi F_ = F.transpose();

	// Send matrix to matlab
	matlab::mlsetmatrix(&engine, "vertex", V_);
	matlab::mlsetmatrix(&engine, "face", F_);
	matlab::mlsetscalar(&engine, "curvature_smoothing", radius);

	//Call matlab engine to calculate
	matlab::mleval(&engine, "[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,curvature_smoothing);");
	
	// Get curvature
	MatrixXd _PD1, _PD2, _PV1, _PV2;
	matlab::mlgetmatrix(&engine, "Umin", _PD1);
	matlab::mlgetmatrix(&engine, "Umax", _PD2);
	matlab::mlgetmatrix(&engine, "Cmin", _PV1);
	matlab::mlgetmatrix(&engine, "Cmax", _PV2);

	PD1 = _PD1.transpose();
	PD2 = _PD2.transpose();
	PV1 = _PV1.col(0);
	PV2 = _PV2.col(0);
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

	Rot = _Rot.block(0, 0, 3, 3);
	//Rot = _Rot;
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

void compute_transformation(const MatrixXi Pair_Index, const MatrixXd V, const VectorXd& PV1, const VectorXd& PV2, const MatrixXd& PD1, const MatrixXd& PD2, const MatrixXd& LN, const int mode, MatrixXd& Transformation)
{
	int size = Pair_Index.rows();
	int dimension;

	switch (mode)
	{
	case TRANSFORM_7D:
		dimension = 7;
		break;
	case TRANSFORM_6D:
		dimension = 6;
		break;
	case ROTATION_3D:
		dimension = 3;
		break;
	case SYMMETRY_PLANE_4D:
		dimension = 4;
		break;
	case ROTATION_AND_SYMMETRY_7D:
		dimension = 7;
		break;
	case SYMMETRY_VECTOR_3D:
		dimension = 3;
		break;
	case ROTATION_AND_SVECTOR_6D:
		dimension = 6;
		break;
	}

	Transformation.resize(size, dimension);
	
	for (int i = 0; i < size; i++)
	{
		MatrixXd Frame_p1, Frame_p2;
		MatrixXd Rot, Rot1, quaternionRot_1, Rot2, quaternionRot_2;
		VectorXd Eular;
		VectorXd Trans;
		VectorXd Plane;
		double Scale;
		Frame_p1.resize(3, 3);
		Frame_p2.resize(3, 3);
		Rot.resize(4, 4);
		Trans.resize(3);
		Plane.resize(4);

		int pair_p1 = Pair_Index(i, 0);
		int pair_p2 = Pair_Index(i, 1);

		if (mode == TRANSFORM_7D)
		{
			Scale = (PV1(pair_p1) / PV1(pair_p2) + PV2(pair_p1) / PV2(pair_p2)) / 2;
		}
		else
		{
			Scale = 1;
		}

		/*Frame_p1.col(0) = LN.row(pair_p1).transpose();
		Frame_p1.col(1) = PD1.row(pair_p1).transpose();
		Frame_p1.col(2) = PD2.row(pair_p1).transpose();

		Frame_p2.col(0) = LN.row(pair_p2).transpose();
		Frame_p2.col(1) = PD1.row(pair_p2).transpose();
		Frame_p2.col(2) = PD2.row(pair_p2).transpose();

		Rot = Frame_p2 * Frame_p1.inverse();
		Trans = V.row(pair_p2).transpose() - Scale * Rot * (V.row(pair_p1).transpose());*/

		//same as paper
		VectorXd p1_p2_n;
		Vector3d origin;
		origin << 0, 0, 0;
		double cos_n1_n2, cos_p1_p2;
		double angle_n1_n2, angle_p1_p2;

		Vector3d N1 = LN.row(pair_p1), N2 = LN.row(pair_p2);
		p1_p2_n = N1.cross(N2);
		p1_p2_n.normalize();
		cos_n1_n2 = (LN.row(pair_p1) * LN.row(pair_p2).transpose())(0,0) / (LN.row(pair_p1).norm()*LN.row(pair_p2).norm());
		angle_n1_n2 = acos(cos_n1_n2);
		compute_rotation(origin, p1_p2_n, angle_n1_n2, Rot1, quaternionRot_1);

		//cout << "rot1" << endl;
		//printMatrix(Rot1);
		//cout << "Q1" << endl;
		//printMatrix(quaternionRot_1);
		//cout << "aim:" << endl;
		//printVector(N2);
		VectorXd N1_ROT;
		N1_ROT = Rot1*N1;
		//cout << "n1 after rot:" << endl;
		//printVector(N1_ROT);



		Vector3d PD1_, PD2_, PD1_ROT;
		PD1_ = (Vector3d)PD1.row(pair_p1);
		PD2_ = (Vector3d)PD1.row(pair_p2);
		PD1_ROT = (Rot1 * PD1_);

		//cout << "PD1 after rot1:" << endl;
		//printVector(PD1_ROT);
		//MatrixXd test;
		//test.resize(3, 3);
		//test = PD1_ROT.transpose()*N1_ROT;
		//cout << "mul N" << endl;
		//printMatrix(test);
		//test = PD1_ROT.transpose()*PD2_;
		//cout << "mul pd2" << endl;
		//printMatrix(test);

		cos_p1_p2 = (PD1_ROT.transpose() * PD2_)(0,0) / (PD1_ROT.norm() * PD2_.norm());
		//bug fixed: cos maybe '>' 1, chair_005 
		cos_p1_p2 = (cos_p1_p2 > 1) ? 1.0:cos_p1_p2;
		//cout << "cos" << cos_p1_p2 << endl;
		angle_p1_p2 = acos(cos_p1_p2);
		VectorXd normal = (PD1_ROT.cross(PD2_));
		//cout << "normal: ";
		//printVector(normal);
		//cout << "N_ROT: ";
		//printVector(N_rot);
		if (normal.dot(N1.head(3)) < 0)
		{
			//cout << "different" << endl;
			angle_p1_p2 = -angle_p1_p2;
		}
		compute_rotation(origin, N1_ROT, angle_p1_p2, Rot2, quaternionRot_2);
		
		//cout << "rot2" << endl;
		//printMatrix(Rot2);
		//cout << "aim:" << endl;
		//printVector(PD2_);
		//Vector3d PD1_ROT_ROT;
		//PD1_ROT_ROT = Rot2*PD1_ROT;
		//cout << "after rot2:" << endl;
		//printVector(PD1_ROT_ROT);
		
		//paper?
		Rot = Rot2 * Rot1;
		//Trans = V.row(pair_p2).transpose() - Scale * Rot * (V.row(pair_p1).transpose());
		//mine
		Trans = V.row(pair_p2).transpose() - Scale * (V.row(pair_p1).transpose());
	
		
		

		
		RotateMatrix2EularAngle(Rot, Eular);
		double A, B, C, D;
		Vector3d midPoint, normalOfPlane;

		midPoint = (V.row(pair_p1) + V.row(pair_p2)) / 2;
		normalOfPlane = V.row(pair_p1) - V.row(pair_p2);
		normalOfPlane.normalize();
		if (normalOfPlane(0) < 0)
		{
			normalOfPlane = normalOfPlane * -1;
		}
		A = normalOfPlane(0);
		B = normalOfPlane(1);
		C = normalOfPlane(2);
		D = -(A*midPoint(0) + B*midPoint(1) + C*midPoint(2));

		switch (mode)
		{
		case TRANSFORM_7D:
			Transformation.row(i) << Scale, Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
			break;
		case TRANSFORM_6D:
			Transformation.row(i) << Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
			break;
		case ROTATION_3D:
			Transformation.row(i) << Eular(0), Eular(1), Eular(2);
			break;
		case SYMMETRY_PLANE_4D:
			Transformation.row(i) << A, B, C, D;
			break;
		case ROTATION_AND_SYMMETRY_7D:
			Transformation.row(i) << Eular(0), Eular(1), Eular(2), A, B, C, D;
			break;
		case SYMMETRY_VECTOR_3D:
			Transformation.row(i) << A, B, C;
			break;
		case ROTATION_AND_SVECTOR_6D:
			Transformation.row(i) << Eular(0), Eular(1), Eular(2), A, B, C;
			break;
		}
		//if (dimension == 7)
		//{
		//	Transformation.row(i) << Scale, Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
		//}
		//else if (dimension == 6)
		//{
		//	Transformation.row(i) << Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
		//}
		//else if (dimension == 4)
		//{
		//	Transformation.row(i) << A, B, C, D;
		//}
		//else if (dimension == 3)
		//{
		//	Transformation.row(i) << Trans(0), Trans(1), Trans(2);
		//}
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
	//if (Rx != Rx)
	//{
	//	cout << "r4, r1 = " << r4 << " , " << r1 << endl;
	//}
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

int MeanShift_Cluster(int mode, const MatrixXd& V, const MatrixXd& Transformation, MatrixXd& clusterCenter, MatrixXd& data2cluster, MatrixXd& cluster2data)
{
	// Send Transformation matrix to matlab
	MatrixXd dataPts_ = Transformation.transpose();
	matlab::mlsetmatrix(&engine, "dataPts", dataPts_);

	matlab::mlsetscalar(&engine, "mode", mode);
	//cout << "mode = " << matlab::mlgetscalar(&engine, "mode") << endl;
	//Cluster using matlab
	matlab::mleval(&engine, "[clustCent,data2cluster,cluster2data,numClust] = MeanShiftCluster(dataPts,mode,boxSize);");
	
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

void patching(const MatrixXd& V, const MatrixXi& F, const MatrixXi& Pair_Index, const MatrixXd& clusterCenter, const MatrixXd& cluster2data, MatrixXd& Patch1, MatrixXd& Patch2, VectorXi& patch_size, VectorXi& patch_rank)
{
	MatrixXd vertex = V.transpose();
	MatrixXi face = F.transpose();
	MatrixXd vertex_in_cluster;

	int cluster_num = clusterCenter.rows();
	int vertex_num = V.rows();

	vertex_in_cluster.setOnes(cluster_num, 2 * vertex_num);
	vertex_in_cluster = -vertex_in_cluster;


	for (int i = 0; i < cluster_num; i++)
	{
		VectorXi vertex_in_this_cluster;
		MatrixXd patch1, patch2;
		Extract_Vertex(cluster2data, Pair_Index, i, vertex_in_this_cluster);
		int vertex_in_this_Cluster_num = vertex_in_this_cluster.rows();
	
		for (int j = 0; j < vertex_in_this_Cluster_num; j++)
		{
			vertex_in_cluster(i, j) = vertex_in_this_cluster(j);
		}
	}
	matlab::mlsetmatrix(&engine, "vertex", vertex);
	matlab::mlsetmatrix(&engine, "face", face);
	matlab::mlsetmatrix(&engine, "transformation", clusterCenter);
	matlab::mlsetmatrix(&engine, "vertex_in_cluster", vertex_in_cluster);
	
	matlab::mleval(&engine, " [ patch1, patch2, patchSize ] = Patching_entrance( vertex, face, transformation, vertex_in_cluster );");

	matlab::mlgetmatrix(&engine, "patch1", Patch1);
	matlab::mlgetmatrix(&engine, "patch2", Patch2);

	MatrixXd patchSize;
	matlab::mlgetmatrix(&engine, "patchSize", patchSize);

	VectorXi temp;
	patch_size.resize(cluster_num);
	for (int j = 0; j < cluster_num; j++)
	{
		patch_size(j) = patchSize(0,j);
	}
	sortrows(patch_size, false, temp, patch_rank);

	//cout << temp;
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

void printVector(const Vector3d& vec)
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

void showCluster(const MatrixXd& Transformation, const MatrixXd& clusterCenter, const MatrixXd& data2cluster)
{
	// Send Transformation matrix to matlab
	MatrixXd dataPts = Transformation.transpose();
	matlab::mlsetmatrix(&engine, "dataPts", dataPts);
	matlab::mlsetmatrix(&engine, "data2cluster", data2cluster);

	//Cluster using matlab
	matlab::mleval(&engine, "[] = showResults(dataPts, data2cluster);");
}
