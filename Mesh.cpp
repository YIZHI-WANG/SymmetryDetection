#include "Mesh.h"

void printMatrix(string name, const MatrixXd& mat)
{
	cout << name << ": " << endl;
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

void printVector(string name, const VectorXd& vec)
{
	cout << name << ": ";
	int size = vec.rows();
	for (int i = 0; i < size; i++)
	{
		cout << vec(i) << ",";
	}
	cout << endl;
}

VectorXd RotateMatrix2EularAngle(const MatrixXd& Rot)
{
	//printMatrix("RotMat", Rot);
	VectorXd Eular;
	Eular.resize(3);
	double
		r1 = Rot(0, 0), r2 = Rot(0, 1), r3 = Rot(0, 2),
		r4 = Rot(1, 0), r5 = Rot(1, 1), r6 = Rot(1, 2),
		r7 = Rot(2, 0), r8 = Rot(2, 1), r9 = Rot(2, 2);
	double Rx = atan2(r4, r1);
	double sx = sin(Rx), cz = cos(Rx);
	double Ry = atan2((-r7), (r1*cz + r4*sx));
	double Rz = atan2((r3*sx - r6*cz), (r5*cz - r2*sx));
	//cout << Rx << "," << Ry << "," << Rz << endl;
	//if (Rx != Rx)
	//{
	//	cout << "r4, r1 = " << r4 << " , " << r1 << endl;
	//}
	Eular << Rx, Ry, Rz;
	//printVector("EularAng", Eular);
	return Eular;
}

MatrixXd compute_rotation(const VectorXd& V_AXIS_start, const VectorXd& V_AXIS_end, const double angle)
{
	//cout << "rot1 angle =" << angle << endl;
	MatrixXd quaternionRot, Rot;
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

	return Rot;
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

Mesh::Mesh()
{
	//init some parameter
	set_sample_ratio(0.5);
	set_prune_threshold(0.75);
	set_pair_threshold(0.1);
	set_smooth_radius(10);
	set_transformation_space_mode(TRANSFORM_6D);
}

Mesh::~Mesh()
{
}

void Mesh::readMesh(string filePath)
{
	string off = "off";
	string obj = "obj";

	transform(filePath.begin(), filePath.end(), filePath.begin(), ::tolower);

	int find_off = filePath.find(off);
	int find_obj = filePath.find(obj);
	if (find_off >= 0)
	{
		readOFF(filePath, vertex, face);
	}
	else if (find_obj >= 0)
	{
		readOBJ(filePath, vertex, face);
	}
	else
	{
		cout << "Sorry, I cannot identify the input file format." << endl;
	}
	vertex_num = vertex.rows();
	face_num = face.rows();
	vertex_num_0 = vertex_num;
	face_num_0 = face_num;
	compute_min_area();

	//sort_vertex_in_face();

	cout << "READ OBJ FINISH!" << endl << endl;
	
	cout << "Model information:" << endl;
	cout << "Vertex: " << vertex_num << endl;
	cout << "Faces: " << face_num << endl;
	cout << endl;
}

void Mesh::set_remesh_fineness(double fineness)
{
	remesh_fineness = fineness;
}

void Mesh::sort_vertex_in_face()
{

	sort_vertex_in_face_by_row();

	sort_vertex_in_face_by_col(0, face_num - 1, 0);


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

void Mesh::sort_vertex_in_face_by_row()
{
	for (int i = 0; i < face_num; i++)
	{
		vector<int> row, row_sorted;
		row.push_back(face(i, 0));
		row.push_back(face(i, 1));
		row.push_back(face(i, 2));
		row_sorted = row;
		sort(row_sorted.begin(), row_sorted.end());
		//cout << "before" << F.row(i) << endl;
		if (row_sorted[0] == row[0])
		{
			face.row(i) << row[0], row[1], row[2];
		}
		else if (row_sorted[0] == row[1])
		{
			face.row(i) << row[1], row[2], row[0];
		}
		else
		{
			face.row(i) << row[2], row[0], row[1];
		}
		//cout << "after" << F.row(i) << endl;
	}

}

void Mesh::sort_vertex_in_face_by_col(int start, int end, int col)
{

	MatrixXi F_ = face;
	VectorXi sortCol, sortedCol, sortedIndex;
	sortCol = F_.block(start, col, end - start + 1, 1);
	sortrows(sortCol, true, sortedCol, sortedIndex);
	for (int i = start; i <= end; i++)
	{
		face.row(i) = F_.row(start + sortedIndex(i - start));
	}
}

void Mesh::remesh()
{
	cout << "Start Remeshing..." << endl;
	sample_aim = ceil(vertex_num_0 * sample_ratio);
	unit_area = totalArea / sample_aim;
	double mean_area = totalArea / 2 / face_num_0 ;
	//double cut_area = (unit_area < min_area_0) ? unit_area : min_area_0;
	//double cut_area = (unit_area < mean_area) ? unit_area : mean_area;
	if (remesh_fineness != 0)
	{
		double cut_area = mean_area / remesh_fineness;
		matlab::mlsetmatrix(&engine, "vertex", vertex);
		matlab::mlsetmatrix(&engine, "face", face);
		matlab::mlsetscalar(&engine, "cut_area", cut_area);

		matlab::mleval(&engine, "[new_vertexs,new_faces] = remesh(vertex, face ,cut_area,1);");

		matlab::mlgetmatrix(&engine, "new_vertexs", vertex);
		matlab::mlgetmatrix(&engine, "new_faces", face);
	}

	vertex_num = vertex.rows();
	face_num = face.rows();

	sort_vertex_in_face();
	compute_area();

	cout << "After Remeshing, " << endl << endl;
	cout << "Vertex: " << vertex_num << endl;
	cout << "Faces: " << face_num << endl;
	cout << endl;
}

void Mesh::compute_local_frame()
{
	compute_boxSize();
	compute_curvature();
	compute_signature();
	compute_local_normal();
}

void Mesh::sample()
{
	cout << "Start Sampling..." << endl;
	
	int sampleNum = 0;
	VectorXi V_Sample;
	V_Sample.resize(vertex_num);
	V_Sample.setZero();


	VectorXd DbGroup_Area;
	DbGroup_Area.resize(vertex_num);
	int currentVertex;
	currentVertex = face(0, 0);
	DbGroup_Area(0) = 0.0;
	for (int i = 0; i < face_num; i++)
	{
		if (face(i, 0) == currentVertex)
		{
			DbGroup_Area(currentVertex) += face_area(i);
		}
		else
		{
			currentVertex = face(i, 0);
			DbGroup_Area(currentVertex) = 0.0;
		}
	}

	//cout << "min = " << min << endl;
	//cout << "max = " << max << endl;

	int sampleUnits = ceil(totalArea / unit_area);
	int samplePerUnit = ceil((double)sample_aim / (double)sampleUnits);
	//cout << sampleUnits << endl;
	//cout << samplePerUnit << endl;
	//cout << "unitArea = " << unitArea << endl;

	double currentGroupArea;
	VectorXi BigFaceIndex;
	int BigFaceNum;
	BigFaceIndex.resize(vertex_num);
	int currentGroupStart, currentGroupEnd;
	currentGroupStart = currentGroupEnd = -1;
	for (int i = 0; i < face_num; i++)
	{
		double DbAreai = face_area(i);
		if (currentGroupStart < 0)
		{
			currentGroupArea = 0;
			currentGroupStart = currentGroupEnd = i;
			BigFaceNum = 0;
			BigFaceIndex.setOnes();
			BigFaceIndex = BigFaceIndex * -1;
		}

		if (DbAreai >= unit_area)
		{
			BigFaceIndex(BigFaceNum++) = i;
			//TODO extract from dbAI
			sample_from_facegroup(samplePerUnit, currentGroupStart, currentGroupEnd, BigFaceNum, BigFaceIndex, sampleNum, V_Sample);

		}
		else
		{
			currentGroupArea += face_area(i);
			currentGroupEnd++;
			if (currentGroupEnd == face_num)
				currentGroupEnd--;
			if (currentGroupArea >= unit_area)
			{
				//cout << "start = " << currentGroupStart << ", end = " << currentGroupEnd << endl;
				//TODO extract from start to end
				sample_from_facegroup(samplePerUnit, currentGroupStart, currentGroupEnd, BigFaceNum, BigFaceIndex, sampleNum, V_Sample);
				currentGroupStart = currentGroupEnd = -1;
			}
		}

	}

	//printVector(V_Sample);

	sample_index.resize(sampleNum);
	int current = 0;
	for (int i = 0; i < vertex_num; i++)
	{
		if (V_Sample(i) == 1)
		{
			sample_index(current++) = i;
		}
	}

	
	cout << "Sampling Finish!" << endl;
	cout << "After Sampling, Vertex: " << sample_index.rows() << ", Sampling Ratio = " << (double)sample_index.rows() / (double)vertex_num << endl << endl;
	prune_umbilic_points();
	cout << "Start Pruning..." << endl;
	cout << "After Pruning, Vertex: " << sample_num << endl << endl;
}

void Mesh::pairing()
{
	cout << "Start Pairing..." << endl;
	sort_sample_by_signature();

	int max_size = sample_num * (sample_num - 1) / 2;
	int _pair_size = 0;

	double dist_threshold = boxSize * 0.03;

	MatrixXi _Pair_Index(max_size, 2);

	if (transformation_space_mode == TRANSFORM_7D || transformation_space_mode == ROT_AND_SYM_AND_SCAL)
	{
		VectorXi floor_index, celling_index;
		find_similar(0,signature.col(0), min_compare2max_sorted, floor_index, celling_index);
		for (int i = 0; i < sample_num; i++)
		{
			if (floor_index(i) != -1 && celling_index(i) != -1)
			{
				//cout << floor_index(i) << "," << celling_index(i) <<endl;
				for (int k = floor_index(i); k <= celling_index(i); k++)
				{
					if (sample_index_dictionary(sample_index_sorted_min_compare2max(k)) != -1 && sample_index_sorted_min_compare2max(k) > pruned_sample_index(i))
					{
						Vector3d dist = vertex.row(pruned_sample_index(i)) - vertex.row(sample_index_sorted_min_compare2max(k));
						if (dist.norm() >= dist_threshold)
						{
							if (signature(pruned_sample_index(i),1) / signature(sample_index_sorted_min_compare2max(k),1) > 0)
							{
								int v1 = pruned_sample_index(i), v2 = sample_index_sorted_min_compare2max(k);
								if (vertex(v1, 0) > vertex(v2, 0) || (vertex(v1, 0) == vertex(v2, 0) && vertex(v1, 1) > vertex(v2, 1)) || (vertex(v1, 0) == vertex(v2, 0) && vertex(v1, 1) == vertex(v2, 1) && vertex(v1, 2) > vertex(v2, 2)))
									//cout << SV(V_Index_Prune(i), 1) << "," << SV(V_Index_Prune(i), 2) << " and "<<SV(V_Index_Sorted_PV1_, 1) << "," << SV(V_Index_Sorted_PV1_, 2) << endl;
								{
									_Pair_Index.row(_pair_size++) << v1, v2;
								}
								else
								{
									_Pair_Index.row(_pair_size++) << v2, v1;
								}
							}
						}
					}
				}
			}
		}
	}
	else
	{
		VectorXi floor_index_PV1, floor_index_PV2, celling_index_PV1, celling_index_PV2;
		find_similar(1,signature.col(1), min_pricinple_val_sorted, floor_index_PV1, celling_index_PV1);
		//printVector("min",min_pricinple_val_sorted);
		//find_similar(V_Index_Prune, SV.col(2), SV_Sorted_PV2, pair_threshold, floor_index_PV2, celling_index_PV2);
		for (int i = 0; i < sample_num; i++)
		{
			int my_Index = pruned_sample_index(i);
			double this_point_pv1 = signature(my_Index, 1);
			double this_point_pv2 = signature(my_Index, 2);
			
			double threshold_value_1 = this_point_pv1 * pair_threshold;
			double threshold_value_2 = this_point_pv1 * pair_threshold;

			double pv1_min = this_point_pv1 - threshold_value_1;
			double pv1_max = this_point_pv1 + threshold_value_1;
			double pv2_min = this_point_pv2 - threshold_value_2;
			double pv2_max = this_point_pv2 + threshold_value_2;


			//if (0< this_point_pv2 && this_point_pv2 <= 0.01)
			//{
			//	pv2_min = 1/INF;
			//	pv2_max = 0.01;
			//}
			//else if (-0.01< this_point_pv2 && this_point_pv2<0 )
			//{
			//	pv2_min = -0.01;
			//	pv2_max = 1/INF;
			//}

			//if (0< this_point_pv2 && this_point_pv2 <= 0.01)
			//{
			//	pv2_min = max(1/INF, this_point_pv2-0.05);
			//	pv2_max = 0.01;
			//}
			//else if (-0.01< this_point_pv2 && this_point_pv2<0 )
			//{
			//	pv2_min = -0.01;
			//	pv2_max = 1/INF;
			//}
				
			
			if (floor_index_PV1(i) != -1 && celling_index_PV1(i) != -1)
			{
				for (int m = floor_index_PV1(i); m <= celling_index_PV1(i); m++)
				{
					int V_Index_Sorted_PV1_ = sample_index_sorted_min_pricin_val(m);
					double pv1 = signature(V_Index_Sorted_PV1_, 1);
					double pv2 = signature(V_Index_Sorted_PV1_, 2);
					if (pv2 >= pv2_min && pv2 <= pv2_max)
					{
						if (sample_index_dictionary(V_Index_Sorted_PV1_) != -1 && V_Index_Sorted_PV1_ > my_Index)
						{
							Vector3d dist = vertex.row(my_Index) - vertex.row(V_Index_Sorted_PV1_);
							if (dist.norm() >= dist_threshold)
							{
								int v1 = my_Index, v2 = V_Index_Sorted_PV1_;
								if (vertex(v1, 0) > vertex(v2, 0) || (vertex(v1, 0) == vertex(v2, 0) && vertex(v1, 1) > vertex(v2, 1)) || (vertex(v1, 0) == vertex(v2, 0) && vertex(v1, 1) == vertex(v2, 1) && vertex(v1, 2) > vertex(v2, 2)))
								{
									_Pair_Index.row(_pair_size++) << v1, v2;
								}
								else
								{
									_Pair_Index.row(_pair_size++) << v2, v1;
								}
								//cout << SV(V_Index_Prune(i), 1) << "," << SV(V_Index_Prune(i), 2) << " and "<<SV(V_Index_Sorted_PV1_, 1) << "," << SV(V_Index_Sorted_PV1_, 2) << endl;
							}
						}
					}
				}
			}
		}
	}
	//cout << "_pair_size = " << _pair_size << endl;
	pairs = _Pair_Index.block(0, 0, _pair_size, 2);
	pair_num = _pair_size;

	cout << "Pairing Finish!" << endl;
	cout << "After Pairing, Pair: " << pair_num << endl << endl;
}

void Mesh::compute_transformation()
{
	cout << "Start Computing Transformation..." << endl;
	int dimension;

	switch (transformation_space_mode)
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
	case ROT_AND_SYM_AND_SCAL:
		dimension = 8;
		break;
	case ROTATION_AND_SVECTOR_6D:
		dimension = 6;
		break;
	}

	transformation.resize(pair_num, dimension);

	for (int i = 0; i < pair_num; i++)
	{
		MatrixXd Frame_p1, Frame_p2;
		MatrixXd Rot, Rot1, Rot2;
		VectorXd Eular;
		VectorXd Trans;
		VectorXd Plane;
		double Scale;
		Frame_p1.resize(3, 3);
		Frame_p2.resize(3, 3);
		Rot.resize(4, 4);
		Trans.resize(3);
		Plane.resize(4);

		int pair_p1 = pairs(i, 0);
		int pair_p2 = pairs(i, 1);

		if (transformation_space_mode == TRANSFORM_7D || transformation_space_mode == ROT_AND_SYM_AND_SCAL)
		{
			Scale = (pricinple_curvanture_min(pair_p1) / pricinple_curvanture_min(pair_p2) + pricinple_curvanture_max(pair_p1) / pricinple_curvanture_max(pair_p2)) / 2;
			//cout << Scale<<":" << pricinple_curvanture_min(pair_p1) <<"," << pricinple_curvanture_min(pair_p2) <<","<< pricinple_curvanture_max(pair_p1) << "," << pricinple_curvanture_max(pair_p2)<< endl;
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

		Vector3d N1 = local_normal.row(pair_p1), N2 = local_normal.row(pair_p2);
		p1_p2_n = N1.cross(N2);
		
		if (p1_p2_n.norm() == 0)
		{
			Rot1.resize(3, 3);
			Rot1.setIdentity();
		/*
			printMatrix("Rot1", Rot1);*/
		}
		else
		{
			p1_p2_n.normalize();
			cos_n1_n2 = (local_normal.row(pair_p1) * local_normal.row(pair_p2).transpose())(0, 0) / (local_normal.row(pair_p1).norm()*local_normal.row(pair_p2).norm());

			cos_n1_n2 = (cos_n1_n2 > 1) ? 1.0 : cos_n1_n2;
			cos_n1_n2 = (cos_n1_n2 < -1) ? -1.0 : cos_n1_n2;
			//cout << "cos1 " << cos_n1_n2 << endl;
			angle_n1_n2 = acos(cos_n1_n2);

			//printVector("p1_p2_n", p1_p2_n);
			Rot1 = compute_rotation(origin, p1_p2_n, angle_n1_n2);
		}
		
		


		VectorXd N1_ROT;
		N1_ROT = Rot1*N1;




		Vector3d PD1_, PD2_, PD1_ROT;
		PD1_ = (Vector3d)pricinple_direction_min.row(pair_p1);
		PD2_ = (Vector3d)pricinple_direction_min.row(pair_p2);
		PD1_ROT = (Rot1 * PD1_);

		//printMatrix("Rot1", Rot1);

		//cout << "PD1_:" << PD1_(0) <<","<< PD1_(1) <<","<< PD1_(2) << endl;


		cos_p1_p2 = (PD1_ROT.transpose() * PD2_)(0, 0) / (PD1_ROT.norm() * PD2_.norm());

		/*cout << "shang " << (PD1_ROT.transpose() * PD2_)(0, 0) << endl;
		cout << "xia" << (PD1_ROT.norm() * PD2_.norm()) << endl;*/
		//bug fixed: cos maybe '>' 1, chair_005 
		cos_p1_p2 = (cos_p1_p2 > 1) ? 1.0 : cos_p1_p2;
		cos_p1_p2 = (cos_p1_p2 < -1) ? -1.0 : cos_p1_p2;
		//cout << "cos2 " << cos_p1_p2 << endl;
		angle_p1_p2 = acos(cos_p1_p2);
		VectorXd normal = (PD1_ROT.cross(PD2_));

		if (normal.dot(N1.head(3)) < 0)
		{
			//cout << "different" << endl;
			angle_p1_p2 = -angle_p1_p2;
		}
		Rot2 = compute_rotation(origin, N1_ROT, angle_p1_p2);


		//paper?
		Rot = Rot2 * Rot1;

		//check reflection: after rotation, pricinple_direction_max are opposite
		bool reflection = false;
		Vector3d PD_max_1, PD_max_2, PD_max_1_Rot;
		PD_max_1 = (Vector3d)pricinple_direction_max.row(pair_p1);
		PD_max_2 = (Vector3d)pricinple_direction_min.row(pair_p2);
		PD_max_1_Rot = (Rot * PD_max_1);
		//if (PD_max_1_Rot.dot(PD_max_2) < 0 && Scale > 0.95 && Scale < 1.05)
		if (PD_max_1_Rot.dot(PD_max_2) < 0)
		{
			reflection = true;
		}

		//Trans = vertex.row(pair_p2).transpose() - Scale * Rot * (vertex.row(pair_p1).transpose());
		//mine
		//Trans = vertex.row(pair_p2).transpose() - Scale * (vertex.row(pair_p1).transpose());
		Trans = vertex.row(pair_p2).transpose() - vertex.row(pair_p1).transpose() * Scale;





		Eular = RotateMatrix2EularAngle(Rot);
		double A, B, C, D;
		Vector3d midPoint, normalOfPlane;

		midPoint = (vertex.row(pair_p1) * Scale  + vertex.row(pair_p2)) / 2;
		normalOfPlane = -Trans.normalized();
		//if (normalOfPlane(0) < 0)
		//{
		//	normalOfPlane = normalOfPlane * -1;
		//	Scale = 1 / Scale;
		//}
		A = normalOfPlane(0);
		B = normalOfPlane(1);
		C = normalOfPlane(2);
		D = -(A*midPoint(0) + B*midPoint(1) + C*midPoint(2));

		switch (transformation_space_mode)
		{
		case TRANSFORM_7D:
			transformation.row(i) << Scale, Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
			break;
		case TRANSFORM_6D:
			transformation.row(i) << Eular(0), Eular(1), Eular(2), Trans(0), Trans(1), Trans(2);
			break;
		case ROTATION_3D:
			transformation.row(i) << Eular(0), Eular(1), Eular(2);
			break;
		case SYMMETRY_PLANE_4D:
			//transformation.row(i) << A, B, C, D;
			if (reflection)
			{
				transformation.row(i) << A, B, C, D;
			}
			else
			{
				transformation.row(i) << 0,0,0,0;
			}
			break;
		case ROTATION_AND_SYMMETRY_7D:
			if (reflection)
			{
				transformation.row(i) << M_PI, M_PI, M_PI, A, B, C, D;
			}
			else
			{
				Vector3d trans_normlized = Trans.normalized();
				transformation.row(i) << Eular(0), Eular(1), Eular(2), trans_normlized(0), trans_normlized(1), trans_normlized(2), Trans.norm();
			}
			
			break;
		case ROT_AND_SYM_AND_SCAL:
			if (reflection)
			{
				transformation.row(i) << Scale, M_PI, M_PI, M_PI, A, B, C, D;
			}
			else
			{
				Vector3d trans_normlized = Trans.normalized();
				transformation.row(i) << Scale, Eular(0), Eular(1), Eular(2), trans_normlized(0), trans_normlized(1), trans_normlized(2), Trans.norm();
			}
			break;
		case ROTATION_AND_SVECTOR_6D:
			transformation.row(i) << Eular(0), Eular(1), Eular(2), A, B, C;
			break;
		}
	}

	cout << "Computing Transformation Finish!" << endl << endl;
}

void Mesh::meanShift_cluster()
{
	cout << "Start Clustering..." << endl;
	// Send Transformation matrix to matlab

	MatrixXd dataPts_ = transformation.transpose();
	int m_mode = (int)transformation_space_mode;
	matlab::mlsetmatrix(&engine, "dataPts", dataPts_);
	matlab::mlsetscalar(&engine, "mode", m_mode);
	matlab::mlsetscalar(&engine, "boxSize", boxSize);
	matlab::mlsetscalar(&engine, "grid", band_width);

	//cout << "mode = " << matlab::mlgetscalar(&engine, "mode") << endl;
	//cout << "mode = " << band_width << endl;
	//Cluster using matlab
	matlab::mleval(&engine, "[clustCent,data2cluster,cluster2data,numClust,clusterSize,~,clusterRank] = MeanShiftCluster(dataPts,mode,boxSize,grid);");

	// Get clustering
	matlab::mlgetmatrix(&engine, "clustCent", cluster_center);
	matlab::mlgetmatrix(&engine, "data2cluster", data2cluster);
	matlab::mlgetmatrix(&engine, "cluster2data", cluster2data);
	matlab::mlgetmatrix(&engine, "clusterSize", cluster_size);
	
	int numClust = (int)matlab::mlgetscalar(&engine, "numClust");
	//double bandWidth = matlab::mlgetscalar(&engine, "bandWidth");
	cluster_num =  numClust;

	MatrixXd rank;
	clusterRank.resize(numClust);
	matlab::mlgetmatrix(&engine, "clusterRank", rank);
	for (int i = 0; i < cluster_num; i++)
	{
		clusterRank(i) = rank(i, 0);
	}

	
	//cout << "bandWidth = " << bandWidth << endl;

	cout << "Clustering Finish!" << endl;
	cout << "After Clustering, Cluster Number: " << cluster_num << endl << endl;
}

void Mesh::extract_vertex_from_cluster(int which_cluster, VectorXi& V_in_Cluster)
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
		V_in_Cluster(i * 2) = pairs(Pair_cluster(i), 0);
		V_in_Cluster(i * 2 + 1) = pairs(Pair_cluster(i), 1);
	}
}

void Mesh::patching()
{
	//cout << "patching_threshold" << patching_threshold << endl;
	cout << "Start Patching..." << endl;
	MatrixXd V = vertex.transpose();
	MatrixXi F = face.transpose();
	MatrixXd vertex_in_cluster;

	vertex_in_cluster.setOnes(cluster_num, 2 * cluster2data.cols());
	vertex_in_cluster = -vertex_in_cluster;


	for (int i = 0; i < cluster_num; i++)
	{
		VectorXi vertex_in_this_cluster;
		MatrixXd patch1, patch2;
		extract_vertex_from_cluster(i, vertex_in_this_cluster);
		int vertex_in_this_Cluster_num = vertex_in_this_cluster.rows();

		for (int j = 0; j < vertex_in_this_Cluster_num; j++)
		{
			vertex_in_cluster(i, j) = vertex_in_this_cluster(j);
		}
	}
	//cout << "patching_threshold" << patching_threshold << endl;
	int patching_mode = (int)transformation_space_mode;
	matlab::mlsetmatrix(&engine, "vertex", V);
	matlab::mlsetmatrix(&engine, "face", F);
	matlab::mlsetmatrix(&engine, "transformation", cluster_center);
	matlab::mlsetmatrix(&engine, "vertex_in_cluster", vertex_in_cluster);
	matlab::mlsetmatrix(&engine, "clusterSize", cluster_size);
	matlab::mlsetscalar(&engine, "mode", patching_mode);
	matlab::mlsetscalar(&engine, "threshold", patching_threshold);
	//cout << "patching_threshold" << patching_threshold << endl;

	//cout << patching_mode << endl;
	matlab::mleval(&engine, " [ patch_count, patch1, patch2, patch_in_cluster_num, patch_2_clsuter,clsuter_2_patch,patchSize] = Patching_entrance( mode, threshold, vertex, face, transformation, vertex_in_cluster, clusterSize);");

	patch_count = matlab::mlgetscalar(&engine, "patch_count");
	matlab::mlgetmatrix(&engine, "patch1", patch1);
	matlab::mlgetmatrix(&engine, "patch2", patch2);

	MatrixXd patch_cluster_num;
	matlab::mlgetmatrix(&engine, "patch_in_cluster_num", patch_cluster_num);
	patch_in_cluster_num.resize(cluster_num);
	for (int i = 0; i < cluster_num; i++)
	{
		int test = patch_cluster_num(0, i);
		patch_in_cluster_num(i) = patch_cluster_num(0, i);
	}

	MatrixXd temp_mat;
	matlab::mlgetmatrix(&engine, "patch_2_clsuter", temp_mat);
	patch_2_clutser.resize(patch_count);
	for (int i = 0; i < patch_count; i++)
	{
		int test = temp_mat(0, i);
		cout << i << "," << test << endl;
 		patch_2_clutser(i) = temp_mat(0, i);
	}
	matlab::mlgetmatrix(&engine, "clsuter_2_patch", clsuter_2_patch);
	
	MatrixXd patchSize;
	matlab::mlgetmatrix(&engine, "patchSize", patchSize);

	VectorXi temp;
	patch_size.resize(patch_count);
	for (int j = 0; j < patch_count; j++)
	{
		patch_size(j) = patchSize(0, j);
	}
	sortrows(patch_size, false, temp, patch_rank);

	//cout << temp;
	cout << "Patching Finish!" << endl << endl;
}

void Mesh::find_similar(int mode,const VectorXd& SV_col, const VectorXd& SV_Sorted, VectorXi& floor_index, VectorXi& celling_index)
{
	int percentage = 2;

	int point_size = SV_col.rows();
	int p_dot_size = sample_num / percentage;
	floor_index.resizeLike(pruned_sample_index);
	celling_index.resizeLike(pruned_sample_index);
	for (int i = 0; i < sample_num; i++)
	{
		floor_index(i) = -1;
		celling_index(i) = -1;
	}

	//TODO: subset of P!
	for (int j = 0; j < p_dot_size; j++)
	{
		int i = j * percentage;
		double _sv = SV_col(pruned_sample_index(i));
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
		/*if (mode == 1)
		{
			if (_sv <= 0.01 && _sv > 0)
			{
				_sv_floor = 1 / INF;
				_sv_celling = 0.01;
			}
			else if (_sv < 0 && _sv >= -0.01)
			{
				_sv_floor = -0.01;
				_sv_celling = -1/INF;
			}
		}*/

		//if (mode == 1)
		//{
		//	if (_sv <= 0.01 && _sv > 0)
		//	{
		//		_sv_floor = max(1 / INF, _sv_floor-0.005);
		//		_sv_celling = _sv_floor+0.05;
		//	}
		//	else if (_sv < 0 && _sv >= -0.01)
		//	{
		//		_sv_floor = _sv_floor - 0.05;
		//		_sv_celling = min(_sv_floor+0.005, -1/INF);
		//	}
		//}

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

void Mesh::sample_from_facegroup(int samplePerUnit, int SampleFaceStart, int SampleFaceEnd, int BigFaceNum, const VectorXi& BigFaceIndex, int& sampleNum, VectorXi& V_Sample)
{
	int thisSampleNum;
	if (SampleFaceStart == SampleFaceEnd)
	{
		/*cout << "here" << endl;*/
		thisSampleNum = (samplePerUnit > 3) ? 3 : samplePerUnit;
		for (int i = 0; i < thisSampleNum; i++)
		{
			int thisPointIndex = face(SampleFaceStart, i);
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
					int thisPointIndex = face(i, j);
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

void Mesh::prune_umbilic_points()
{
	VectorXi _V_Index_Prune;
	int size_Sample = sample_index.rows();
	int size_P = 0;
	_V_Index_Prune.resize(size_Sample, 1);
	sample_index_dictionary.resize(vertex_num, 1);
	sample_index_dictionary.fill(-1);

	//cout << "size = " << size << endl;
	for (int i = 0; i < size_Sample; i++)
	{
		//cout << i << "i = " << SV(i) << endl;
		//cout << fabs(signature(sample_index(i), 0) << endl;
		//cout << prune_threshold << endl;
		if (signature(sample_index(i), 1) >0 && signature(sample_index(i), 2) > 0 && signature(sample_index(i),0) < prune_threshold)
		{
			/*cout << "1 " <<V_Index_Sample(i) << endl;
			cout << "2 " <<V_Index_Dic(V_Index_Sample(i)) << endl;*/
			sample_index_dictionary(sample_index(i)) = size_P;
			_V_Index_Prune(size_P++) = sample_index(i);
			//cout << "add " << i << endl;
		}
		else if(signature(sample_index(i), 1) / (signature(sample_index(i), 2)+0.000000001) < 0)
		{
			sample_index_dictionary(sample_index(i)) = size_P;
			_V_Index_Prune(size_P++) = sample_index(i);
		}
		else if (signature(sample_index(i), 1) <0 && signature(sample_index(i), 2) < 0 && signature(sample_index(i), 0) > 1/prune_threshold)
		{
			sample_index_dictionary(sample_index(i)) = size_P;
			_V_Index_Prune(size_P++) = sample_index(i);
		}
		//else
		//{
		//	cout << signature(sample_index(i), 1) <<","<< signature(sample_index(i), 2) << endl;
		//}
	}

	pruned_sample_index = _V_Index_Prune.head(size_P);
	sample_num = size_P;
}

void Mesh::sort_sample_by_signature()
{
	VectorXd SV_7d = signature.col(0);
	VectorXd SV_PV1 = signature.col(1);
	VectorXd SV_PV2 = signature.col(2);
	sortrows(SV_7d, true, min_compare2max_sorted, sample_index_sorted_min_compare2max);
	sortrows(SV_PV1, true, min_pricinple_val_sorted, sample_index_sorted_min_pricin_val);
	sortrows(SV_PV2, true, max_pricinple_val_sorted, sample_index_sorted_max_pricin_val);
}

void Mesh::compute_boxSize()
{
	// Send V matrix to matlab for bounding box size
	MatrixXd V_ = vertex.transpose();
	matlab::mlsetmatrix(&engine, "V", V_);
	matlab::mleval(&engine, "[boxSize] = BoundingBoxSize(V);");

	boxSize = matlab::mlgetscalar(&engine, "boxSize");

	cout << "boxSize = " << boxSize << endl << endl;
}

void Mesh::compute_curvature()
{
	MatrixXd V_ = vertex.transpose();
	MatrixXi F_ = face.transpose();
	double radius = smooth_radius;
	// Send matrix to matlab
	matlab::mlsetmatrix(&engine, "vertex", V_);
	matlab::mlsetmatrix(&engine, "face", F_);
	matlab::mlsetscalar(&engine, "curvature_smoothing", radius);
	matlab::mleval(&engine, "options.curvature_smoothing = curvature_smoothing");

	//Call matlab engine to calculate
	matlab::mleval(&engine, "[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,face,options);");

	// Get curvature
	MatrixXd _PD1, _PD2, _PV1, _PV2;
	matlab::mlgetmatrix(&engine, "Umin", _PD1);
	matlab::mlgetmatrix(&engine, "Umax", _PD2);
	matlab::mlgetmatrix(&engine, "Cmin", _PV1);
	matlab::mlgetmatrix(&engine, "Cmax", _PV2);

	pricinple_direction_min = _PD1.transpose();
	pricinple_direction_max = _PD2.transpose();
	pricinple_curvanture_min = _PV1.col(0);
	pricinple_curvanture_max = _PV2.col(0);
}

void Mesh::compute_signature()
{
	signature.resize(vertex_num, 3);

	for (int i = 0; i < vertex_num; i++)
	{
		if (pricinple_curvanture_max(i) == 0.0)
		{
			signature(i, 0) = (pricinple_curvanture_min(i) > 0.0) ? (-INF) : (INF);
		}
		else
		{
			signature(i, 0) = pricinple_curvanture_min(i) / pricinple_curvanture_max(i);
		}
		//cout << "SV" << i << " = " << SV(i) << endl;
	}

	signature.col(1) = pricinple_curvanture_min;
	signature.col(2) = pricinple_curvanture_max;
}

void Mesh::compute_local_normal()
{
	local_normal.resize(vertex_num, 3);
	for (int i = 0; i < vertex_num; i++)
	{
		Vector3d d1 = pricinple_direction_min.row(i);
		Vector3d d2 = pricinple_direction_max.row(i);
		local_normal.row(i) = d1.cross(d2);
	}
}

void Mesh::compute_area()
{
	doublearea(vertex, face, face_area);
	totalArea = 0.0;
	for (int i = 0; i < face_num; i++)
	{
		totalArea += face_area(i);
	}
}

void Mesh::compute_min_area()
{
	doublearea(vertex, face, face_area);
	min_area_0 = INF, max_area_0 = 0;
	for (int i = 0; i < face_num; i++)
	{
		if (min_area_0 > face_area(i))
			min_area_0 = face_area(i);
		if (max_area_0 < face_area(i))
			max_area_0 = face_area(i);
	}
}

void Mesh::set_smooth_radius(double radius)
{
	smooth_radius = radius;
}

void Mesh::set_sample_ratio(double ratio)
{
	sample_ratio = ratio;
}

void Mesh::set_prune_threshold(double threshold)
{
	prune_threshold = threshold;
}

void Mesh::set_pair_threshold(double threshold)
{
	pair_threshold = threshold;
}

void Mesh::set_band_width(double bandwidth)
{
	band_width = bandwidth;
}

void Mesh::set_patching_threshold(double threshold)
{
	patching_threshold = threshold;
}

void Mesh::set_transformation_space_mode(mode value)
{
	transformation_space_mode = value;
}

void Mesh::init_data()
{
	pricinple_direction_min.resize(0, 0);
	pricinple_direction_max.resize(0, 0);
	local_normal.resize(0, 0);
	signature.resize(0, 0);
	pairs.resize(0, 0);
	transformation.resize(0, 0);
	cluster_center.resize(0, 0);
	data2cluster.resize(0, 0);
	cluster2data.resize(0, 0);
	patch1.resize(0, 0);
	patch2.resize(0, 0);
	signature.resize(0, 0);
	pricinple_curvanture_min.resize(0);
	pricinple_curvanture_max.resize(0);
	face_area.resize(0);
	sample_index.resize(0);
	pruned_sample_index.resize(0);
	sample_index_dictionary.resize(0);
	min_compare2max_sorted.resize(0);
	min_pricinple_val_sorted.resize(0);
	max_pricinple_val_sorted.resize(0);
	sample_index_sorted_min_compare2max.resize(0);
	sample_index_sorted_min_pricin_val.resize(0);
	sample_index_sorted_max_pricin_val.resize(0);
	patch_size.resize(0);
	patch_rank.resize(0);
	compute_area();
	start_engine();
}

void Mesh::start_engine()
{
	cout << "Start Lunching Matlab Engine..." << endl;
	// Launch MATLAB
	matlab::mlinit(&engine);
	if (engine == 0)
	{
		cout << "Lunch Matlab Failed!" << endl;
	}
	else
	{
		cout << "Lunch Matlab Success!" << endl << endl;
		matlab::mleval(&engine, "clear all");
	}
}

bool Mesh::close_engine()
{
	// Close MATLAB
	matlab::mlclose(&engine);
	return true;
}

const MatrixXd& Mesh::get_vertex()
{
	return vertex;
}

const MatrixXi& Mesh::get_face()
{
	return face;
}

const MatrixXd& Mesh::get_signature()
{
	return signature;
}

const VectorXi & Mesh::get_sample_index()
{
	return sample_index;
}

const VectorXi & Mesh::get_pruned_index()
{
	return pruned_sample_index;
}

const MatrixXi & Mesh::get_pairs()
{
	return pairs;
}

VectorXd Mesh::get_cluster_center(int which_cluster)
{
	int dimension = cluster_center.cols();
	VectorXd center;
	center.resize(dimension);
	center = cluster_center.row(which_cluster);
	return center;
}

VectorXd Mesh::get_patch(int which, int which_cluster, int which_patch)
{
	int patch_row = clsuter_2_patch(which_cluster, which_patch);
	VectorXd patch;
	if (patch_row >= 0)
	{
		if (which == 1)
		{

			patch = patch1.row(patch_row);
		}
		else
		{
			patch = patch2.row(patch_row);
		}
		return patch;
	}
	else
	{
		patch.resize(0);
		return patch;
	}
}

VectorXd Mesh::get_patch(int which, int patch_row)
{
	VectorXd patch;
	if (patch_row >= 0)
	{
		if (which == 1)
		{

			patch = patch1.row(patch_row);
		}
		else
		{
			patch = patch2.row(patch_row);
		}
		return patch;
	}
	else
	{
		patch.resize(0);
		return patch;
	}
}
int Mesh::get_nth_patch_index(int n)
{
	return patch_rank(n);
}

int Mesh::get_patch_size(int n)
{
	return patch_size(n);
}

int Mesh::get_cluster_num()
{
	return cluster_num;
}

int Mesh::get_max_patch(int cluster)
{
	return patch_in_cluster_num(cluster);
}

int Mesh::get_cluster_from_rank(int rank)
{
	return clusterRank(rank);
}

int Mesh::get_cluster_from_patch(int which_patch)
{
	return patch_2_clutser(which_patch);
}

int Mesh::get_patch_num()
{
	return patch_count;
}

double Mesh::get_boxSize()
{
	return boxSize;
}

