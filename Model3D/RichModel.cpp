#include "RichModel.h"
#include "EdgePoint.h"
#include <queue>
#include <cassert>
#include <math.h>
#include <iostream>
#include <queue>
#include <functional>
#define _USE_MATH_DEFINES
#include "Parameters.h"
#include "CDT.h"
using namespace std;
namespace Model3D
{
	CRichModel::CRichModel(const CBaseModel &model) : CBaseModel(model)
	{
		PreprocessBaseModelIntoRichModel();
	}

	CRichModel::CRichModel(const string& filename) : CBaseModel(filename)
	{
	}

	CRichModel::CRichModel(const vector<CPoint3D> &verts, const vector<CBaseModel::CFace> &faces) : CBaseModel("")
	{
		m_Verts = verts;
		m_Faces = faces;
		CBaseModel::PreprocessVertsAndFacesIntoBaseModel();
		PreprocessBaseModelIntoRichModel();
	}

	void CRichModel::CreateEdgesFromVertsAndFaces()
	{
		m_Edges.clear();
		m_Edges.reserve(2 * (GetNumOfVerts() + GetNumOfFaces() - 2));
		map<pair<int, int>, int> pondOfUndeterminedEdges;
		int szFaces = GetNumOfFaces();
		for (int i = 0; i < szFaces; ++i)
		{
			int threeIndices[3];
			for (int j = 0; j < 3; ++j)
			{
				int post = (j + 1) % 3;
				int pre = (j + 2) % 3;

				int leftVert = Face(i)[pre];
				int rightVert = Face(i)[j];

				map<pair<int, int>, int>::const_iterator it = pondOfUndeterminedEdges.find(make_pair(leftVert, rightVert));
				if (it != pondOfUndeterminedEdges.end())
				{
					int posInEdgeList = it->second;
					if (m_Edges[posInEdgeList].indexOfOppositeVert != -1)
					{
						ofstream out("tmp.obj");
						out << "g 3dLine\n";
						out << "v " << Vert(m_Edges[posInEdgeList].indexOfLeftVert).x
							<< " " << Vert(m_Edges[posInEdgeList].indexOfLeftVert).y
							<< " " << Vert(m_Edges[posInEdgeList].indexOfLeftVert).z << endl;
						out << "v " << Vert(m_Edges[posInEdgeList].indexOfRightVert).x
							<< " " << Vert(m_Edges[posInEdgeList].indexOfRightVert).y
							<< " " << Vert(m_Edges[posInEdgeList].indexOfRightVert).z << endl;
						out << "l 1 2\n";
						out.close();
						cerr << "Repeated edges!" << endl;
						//cerr << Vert(m_Edges[posInEdgeList].indexOfLeftVert) << endl;
						throw "Repeated edges!";
					}
					threeIndices[j] = posInEdgeList;
					m_Edges[posInEdgeList].indexOfOppositeVert = Face(i)[post];
					m_Edges[posInEdgeList].indexOfFrontFace = i;
				}
				else
				{
					CEdge edge;
					edge.indexOfLeftVert = leftVert;
					edge.indexOfRightVert = rightVert;
					edge.indexOfFrontFace = i;
					edge.indexOfOppositeVert = Face(i)[post];
					edge.indexOfReverseEdge = (int)m_Edges.size() + 1;
					edge.length = (Vert(leftVert) - Vert(rightVert)).Len();
					m_Edges.push_back(edge);
					pondOfUndeterminedEdges[make_pair(leftVert, rightVert)] = threeIndices[j] = (int)m_Edges.size() - 1;

					edge.indexOfLeftVert = rightVert;
					edge.indexOfRightVert = leftVert;
					edge.indexOfReverseEdge = (int)m_Edges.size() - 1;
					edge.indexOfOppositeVert = -1;
					edge.indexOfFrontFace = -1;
					m_Edges.push_back(edge);
					pondOfUndeterminedEdges[make_pair(rightVert, leftVert)] = (int)m_Edges.size() - 1;
				}
			}
			for (int j = 0; j < 3; ++j)
			{
				m_Edges[threeIndices[j]].indexOfLeftEdge = Edge(threeIndices[(j + 2) % 3]).indexOfReverseEdge;
				m_Edges[threeIndices[j]].indexOfRightEdge = Edge(threeIndices[(j + 1) % 3]).indexOfReverseEdge;
			}
		}
		m_Edges.swap(vector<CEdge>(m_Edges));
	}

	void CRichModel::CollectAndArrangeNeighs()
	{
		m_nIsolatedVerts = 0;
		vector<int> sequenceOfDegrees(GetNumOfVerts(), 0);
		m_NeighsAndAngles.clear();
		m_NeighsAndAngles.resize(GetNumOfVerts());
		for (int i = 0; i < (int)m_NeighsAndAngles.size(); ++i)
		{
			m_NeighsAndAngles[i].resize(1, make_pair(-1, 0));
		}
		for (int i = 0; i < (int)GetNumOfEdges(); ++i)
		{
			const CEdge& edge = Edge(i);
			++sequenceOfDegrees[edge.indexOfLeftVert];
			int &indexOfStartEdge = m_NeighsAndAngles[edge.indexOfLeftVert][0].first;
			if (indexOfStartEdge == -1 || !IsStartEdge(indexOfStartEdge))
			{
				indexOfStartEdge = i;
			}
			else if (IsStartEdge(i))
			{
				m_NeighsAndAngles[edge.indexOfLeftVert].push_back(make_pair(i, 0));
			}
		}
		for (int i = 0; i < GetNumOfVerts(); ++i)
		{
			if (m_NeighsAndAngles[i][0].first == -1)
			{
				m_NeighsAndAngles[i].clear();
				m_nIsolatedVerts++;
				continue;
			}
			vector<int> startEdges;
			for (int j = 0; j < (int)Neigh(i).size(); ++j)
			{
				startEdges.push_back(Neigh(i)[j].first);
			}
			m_NeighsAndAngles[i].resize(sequenceOfDegrees[i], make_pair(0, 0));
			int num(0);
			for (int j = 0; j < (int)startEdges.size(); ++j)
			{
				int curEdge = startEdges[j];
				while (1)
				{
					m_NeighsAndAngles[i][num].first = curEdge;
					++num;
					if (num >= sequenceOfDegrees[i])
						break;
					if (IsExtremeEdge(curEdge))
						break;
					curEdge = Edge(curEdge).indexOfLeftEdge;
					if (curEdge == startEdges[j])
					{
						break;
					}
				}
			}
			if (num != sequenceOfDegrees[i])
			{
				throw "Complex vertices";
			}
		}
	}

	void CRichModel::ComputeAnglesAroundVerts()
	{
		m_FlagsForCheckingConvexVerts.clear();
		m_FlagsForCheckingConvexVerts.resize(GetNumOfVerts());
		//for (int i = 0; i < (int)m_NeighsAndAngles.size(); ++i)
		//{
		//	m_NeighsAndAngles[i].resize(Neigh(i).size());
		//}
		for (int i = 0; i < (int)m_NeighsAndAngles.size(); ++i)
		{
			double angleSum(0);
			for (int j = 0; j < (int)m_NeighsAndAngles[i].size(); ++j)
			{
				if (IsExtremeEdge(Neigh(i)[j].first))
					m_NeighsAndAngles[i][j].second = 2 * M_PI + 0.1;
				else
				{
					int next = j + 1;
					if (next >= (int)m_NeighsAndAngles[i].size())
					{
						next = 0;
					}
					double l = Edge(Neigh(i)[j].first).length;
					double r = Edge(Neigh(i)[next].first).length;
					double b = Edge(Edge(Neigh(i)[j].first).indexOfRightEdge).length;
					m_NeighsAndAngles[i][j].second = acos((l * l + r * r - b * b) / (2 * l * r));
					m_Edges[Edge(Edge(Neigh(i)[j].first).indexOfRightEdge).indexOfReverseEdge].angleOpposite = m_NeighsAndAngles[i][j].second;
				}
				angleSum += m_NeighsAndAngles[i][j].second;
			}
			m_FlagsForCheckingConvexVerts[i].first = (angleSum < 2 * M_PI - 5 * AngleTolerance);
			m_FlagsForCheckingConvexVerts[i].second = (angleSum < 2 * M_PI + 5 * AngleTolerance);
		}

	}

	void CRichModel::ComputePlanarCoordsOfIncidentVertForEdges()
	{
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (IsExtremeEdge(i))
				continue;
			double bottom = Edge(i).length;
			double leftLen = Edge(Edge(i).indexOfLeftEdge).length;
			double squareOfLeftLen = leftLen * leftLen;
			double rightLen = Edge(Edge(i).indexOfRightEdge).length;
			double x = (squareOfLeftLen - rightLen * rightLen) / bottom + bottom;
			x /= 2.0;
			m_Edges[i].coordOfOppositeVert = Eigen::Vector2d(x, sqrt(max(0.0, squareOfLeftLen - x * x)));
		}
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (IsExtremeEdge(i))
				continue;
			{
				int reverseEdge = m_Edges[m_Edges[i].indexOfLeftEdge].indexOfReverseEdge;
				Eigen::Vector2d coord = GetNew2DCoordinatesByReversingCurrentEdge(reverseEdge, m_Edges[reverseEdge].coordOfOppositeVert);
				double scale = abs(coord(0)) + abs(coord(1));
				coord(0) /= scale;
				coord(1) /= scale;
				double len = sqrt(coord(0) * coord(0) + coord(1) * coord(1));
				m_Edges[i].matrixRotatedToLeftEdge = Eigen::Vector2d(coord(0) / len, coord(1) / len);
			}
			{
				int reverseEdge = m_Edges[m_Edges[i].indexOfRightEdge].indexOfReverseEdge;
				double rightX = m_Edges[reverseEdge].length;
				double rightY = 0;
				double leftX = m_Edges[reverseEdge].length - m_Edges[reverseEdge].coordOfOppositeVert(0);
				double leftY = -m_Edges[reverseEdge].coordOfOppositeVert(1);

				double detaX = rightX - leftX;
				double detaY = rightY - leftY;
				double scale = abs(detaX) + abs(detaY);
				detaX /= scale;
				detaY /= scale;
				double len = sqrt(detaX * detaX + detaY * detaY);
				m_Edges[i].matrixRotatedToRightEdge = Eigen::Vector2d(detaX / len, detaY / len);
			}
		}
	}

	void CRichModel::FinishChangingEdgeLengths()
	{
		m_maxEdgeLength = 0;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (Edge(i).length > m_maxEdgeLength)
				m_maxEdgeLength = Edge(i).length;
		}

		//Perform the following 
		//when edge lengths are changed.
		//compute angles around a vertex
		ComputeAnglesAroundVerts();
		//planar unfolding
		ComputePlanarCoordsOfIncidentVertForEdges();
	}

	//void CRichModel::MakeOpen2Closed()
	//{
	//	set<int> extremeEdges;
	//	for (int i = 0; i < GetNumOfEdges(); ++i)
	//	{
	//		if (IsExtremeEdge(i))
	//			extremeEdges.insert(i);
	//	}
	//	while (!extremeEdges.empty())
	//	{
	//		vector<int> vertSequence;
	//		int edge = *extremeEdges.begin();
	//		do
	//		{
	//			int v = Edge(edge).indexOfLeftVert;
	//			vertSequence.push_back(v);
	//			extremeEdges.erase(edge);
	//			int nxt = Edge(edge).indexOfRightVert;
	//			edge = Neigh(nxt).back().first;
	//			if (extremeEdges.find(edge) == extremeEdges.end())
	//				break;
	//		} while (true);
	//		for (int i = 1; i < vertSequence.size()-1; ++i)
	//		{
	//			m_Faces.push_back(CFace(vertSequence[0], vertSequence[i], vertSequence[i+1]));
	//		}
	//	}
	//	CBaseModel::PreprocessVertsAndFacesIntoBaseModel();
	//	PreprocessBaseModelIntoRichModel();
	//}

	map<CPoint3D, int> CRichModel::AddFaceInteriorPoints(const vector<pair<int, CPoint3D>> &pts)
	{
		auto res = CBaseModel::AddFaceInteriorPoints(pts);
		PreprocessBaseModelIntoRichModel();
		return res;
	}

	int CRichModel::AddFaceInteriorPoints(int faceIndex, CPoint3D pt)
	{
		int id = CBaseModel::AddFaceInteriorPoints(faceIndex, pt);
		PreprocessBaseModelIntoRichModel();
		return id;
	}

	void CRichModel::MakeOpen2Closed()
	{
		set<int> extremeEdges;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (IsExtremeEdge(i))
				extremeEdges.insert(i);
		}
		while (!extremeEdges.empty())
		{
			vector<int> vertSequence;
			int edge = *extremeEdges.begin();
			do
			{
				int v = Edge(edge).indexOfLeftVert;
				vertSequence.push_back(v);
				extremeEdges.erase(edge);
				int nxt = Edge(edge).indexOfRightVert;
				edge = Neigh(nxt).back().first;
				if (extremeEdges.find(edge) == extremeEdges.end())
					break;
			} while (true);

			CPoint3D normal(0, 0, 0);
			for (int i = 1; i < vertSequence.size() - 1; ++i)
			{
				normal = normal + VectorCross(Vert(vertSequence[0]), Vert(vertSequence[i]), Vert(vertSequence[i + 1]));
			}
			normal.Normalize();
			auto dir1 = normal.GetUnitPerpendicularDir();
			auto dir2 = normal * dir1;
			dir2.Normalize();
			vector<Point_2> pts;
			map<Point_2, int> pts2ID;
			for (int i = 0; i < vertSequence.size(); ++i)
			{
				auto pt3d = Vert(vertSequence[i]);
				Point_2 pt2d(pt3d ^ dir1, pt3d ^ dir2);
				pts.push_back(pt2d);
				pts2ID[pt2d] = vertSequence[i];
			}
			CDT cdt;
			cdt.insert_constraint(pts.begin(), pts.end(), true);
			mark_domains(cdt);
			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
				fit != cdt.finite_faces_end(); ++fit)
			{
				if (fit->info().in_domain())
				{
					int id1 = pts2ID[fit->vertex(0)->point()];
					int id2 = pts2ID[fit->vertex(1)->point()];
					int id3 = pts2ID[fit->vertex(2)->point()];
					//cerr << "ids : " << id1 << ", " << id2 << ", " << id3 << endl;
					m_Faces.push_back(CFace(id1, id2, id3));
				}
			}
		}
		CBaseModel::PreprocessVertsAndFacesIntoBaseModel();
		PreprocessBaseModelIntoRichModel();
	}

	//for (CDT::Finite_faces_iterator fit = denseCDT.finite_faces_begin();
	//	fit != denseCDT.finite_faces_end(); ++fit)
	//{
	//	if (fit->info().in_domain())
	//	{
	//		Polygon_2 p;
	//		p.push_back(fit->vertex(0)->point());
	//		p.push_back(fit->vertex(1)->point());
	//		p.push_back(fit->vertex(2)->point());
	//		linearApprox[fit] = MA::Linear_function<K>(fit->vertex(0)->point(), f(fit->vertex(0)->point()),
	//			fit->vertex(1)->point(), f(fit->vertex(1)->point()),
	//			fit->vertex(2)->point(), f(fit->vertex(2)->point()));
	//	}
	//}

	void CRichModel::PreprocessBaseModelIntoRichModel()
	{
		m_coord_systems_faces.clear();
		for (int i = 0; i < GetNumOfFaces(); ++i)
		{
			auto dir1 = Vert(Face(i)[1]) - Vert(Face(i)[0]);
			dir1.Normalize();
			auto normal = VectorCross(Vert(Face(i)[0]), Vert(Face(i)[1]), Vert(Face(i)[2]));
			normal.Normalize();
			auto y = normal * dir1;
			y.Normalize();
			m_coord_systems_faces.push_back(make_tuple(Vert(Face(i)[0]), dir1, y, normal));
		}
		m_UselessEdges.clear();
		//build edges and compute lengths
		//cerr << "Line " << __LINE__ << endl;
		CreateEdgesFromVertsAndFaces();
		m_maxEdgeLength = 0;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (Edge(i).length > m_maxEdgeLength)
				m_maxEdgeLength = Edge(i).length;
		}
		//cerr << "Line " << __LINE__ << endl;
		//build edges incident to a vertex
		CollectAndArrangeNeighs();
		//cerr << "Line " << __LINE__ << endl;
		//num of open boundaries
		ComputeNumOfHoles();
		//num of components
		//cerr << "Line " << __LINE__ << endl;
		ComputeNumOfComponents();
		//cerr << "Line " << __LINE__ << endl;
		//Perform the following 
		//when edge lengths are changed.
		//compute angles around a vertex
		ComputeAnglesAroundVerts();
		//planar unfolding
		//cerr << "Line " << __LINE__ << endl;
		ComputePlanarCoordsOfIncidentVertForEdges();

		ExtractBoundary();
	}

	void CRichModel::LoadModel()
	{
		CBaseModel::LoadModel();
		PreprocessBaseModelIntoRichModel();
	}


	void CRichModel::ComputeNumOfHoles()
	{
		m_nBoundries = 0;
		if (IsClosedModel())
		{
			return;
		}
		set<int> extremeEdges;
		for (int i = 0; i < (int)m_Edges.size(); ++i)
		{
			if (m_Edges[i].indexOfOppositeVert != -1)
				continue;
			extremeEdges.insert(i);
		}

		while (!extremeEdges.empty())
		{
			++m_nBoundries;
			int firstEdge = *extremeEdges.begin();
			int edge = firstEdge;
			do
			{
				extremeEdges.erase(edge);
				int root = Edge(edge).indexOfRightVert;
				int index = GetSubindexToVert(root, Edge(edge).indexOfLeftVert);
				edge = Neigh(root)[(index - 1 + (int)Neigh(root).size()) % (int)Neigh(root).size()].first;
			} while (edge != firstEdge && !extremeEdges.empty());
		}
	}

	void CRichModel::ComputeNumOfComponents()
	{
		m_nComponents = 0;
		vector<bool> flags(GetNumOfVerts(), false);
		int cnt(0);
		while (cnt < GetNumOfVerts())
		{
			int v;
			for (int i = 0; i < (int)flags.size(); ++i)
			{
				if (!flags[i])
				{
					v = i;
					break;
				}
			}
			queue<int> Que;
			Que.push(v);
			while (!Que.empty())
			{
				int v = Que.front();
				Que.pop();
				if (flags[v])
					continue;
				flags[v] = true;
				cnt++;
				for (int i = 0; i < (int)Neigh(v).size(); ++i)
				{
					if (!flags[Edge(Neigh(v)[i].first).indexOfRightVert])
					{
						Que.push(Edge(Neigh(v)[i].first).indexOfRightVert);
					}
				}
			}
			m_nComponents++;
		}
	}

	void CRichModel::PrintInfo(ostream& out) const
	{
		out << "Model info is as follows.\n";
		out << "Name: " << GetFileShortName() << endl;
		out << "VertNum = " << GetNumOfVerts() << endl;
		out << "FaceNum = " << GetNumOfFaces() << endl;
		out << "EdgeNum = " << GetNumOfEdges() << endl;
		out << "Scale = " << m_scale << endl;
		if (!IsClosedModel())
			out << "BoundaryNum = " << GetNumOfBoundries() << endl;
		out << "Genus = " << GetNumOfGenera() << endl;
		if (IsClosedModel())
			out << "It is a closed model.\n";
		else
			out << "It is an open model.\n";
		if (GetNumOfComponents() != 1)
			out << "It has " << GetNumOfComponents() << " components.\n";
	}

	void CRichModel::SetEdgeLength(int leftVert, int rightVert, double newLength)
	{
		int edgeID = GetEdgeIndexFromTwoVertices(leftVert, rightVert);
		int reverseID = Edge(edgeID).indexOfReverseEdge;
		m_Edges[edgeID].length = newLength;
		m_Edges[reverseID].length = newLength;
	}


	double CRichModel::AngleSum(int vertIndex) const
	{
		double angleSum(0);
		for (int j = 0; j < (int)m_NeighsAndAngles[vertIndex].size(); ++j)
		{
			angleSum += m_NeighsAndAngles[vertIndex][j].second;
		}
		return angleSum;
	}

	pair<double, double> CRichModel::GetTwoSplitAngles(int root, EdgePoint pt1, EdgePoint pt2) const
	{
		if (!pt1.isVertex && pt2.isVertex)
		{
			pair<double, double> splitAngles = GetTwoSplitAngles(root, pt2, pt1);
			return make_pair(splitAngles.second, splitAngles.first);
		}
		if (!pt1.isVertex)
		{
			if (Edge(pt1.index).indexOfOppositeVert != root)
			{
				pt1.proportion = 1 - pt1.proportion;
				pt1.index = Edge(pt1.index).indexOfReverseEdge;
				assert(Edge(pt1.index).indexOfOppositeVert == root);
			}
		}
		if (!pt2.isVertex)
		{
			if (Edge(pt2.index).indexOfOppositeVert != root)
			{
				pt2.proportion = 1 - pt2.proportion;
				pt2.index = Edge(pt2.index).indexOfReverseEdge;
				assert(Edge(pt2.index).indexOfOppositeVert == root);
			}
		}

		double rightAngleSum(0);
		double leftAngleSum(0);
		if (pt1.isVertex && pt2.isVertex)
		{
			int index1 = GetSubindexToVert(root, pt1.index);
			int index2 = GetSubindexToVert(root, pt2.index);
			int index = index1;
			while (index != index2)
			{
				rightAngleSum += Neigh(root)[index].second;
				index = (index + 1) % Neigh(root).size();
			}
			index = index2;
			leftAngleSum = Neigh(root)[index].second;
			index = (index + 1) % Neigh(root).size();
			while (index != index1)
			{
				leftAngleSum += Neigh(root)[index].second;
				index = (index + 1) % Neigh(root).size();
			}
		}
		else if (pt1.isVertex && !pt2.isVertex)
		{
			int index1 = GetSubindexToVert(root, pt1.index);
			int index2 = GetSubindexToVert(root, Edge(pt2.index).indexOfLeftVert);
			int index = index1;
			while (index != index2)
			{
				rightAngleSum += Neigh(root)[index].second;
				index = (index + 1) % Neigh(root).size();
			}
			int index3 = GetSubindexToVert(root, Edge(pt2.index).indexOfRightVert);
			index = index3;
			while (index != index1)
			{
				leftAngleSum += Neigh(root)[index].second;
				index = (index + 1) % Neigh(root).size();
			}

			if (pt2.proportion < 1e-2
				|| pt2.proportion > 1 - 1e-2
				|| Edge(pt2.index).angleOpposite < 5.0 * M_PI / 180)
			{
				double leftAngle2 = pt2.proportion * Edge(pt2.index).angleOpposite;
				double rightAngle2 = (1 - pt2.proportion) * Edge(pt2.index).angleOpposite;
				rightAngleSum += leftAngle2;
				leftAngleSum += rightAngle2;
			}
			else
			{
				double c = Edge(pt2.index).length * pt2.proportion;
				double a = Edge(Edge(pt2.index).indexOfLeftEdge).length;
				double detaX = c - Edge(pt2.index).coordOfOppositeVert(0);
				double detaY = 0 - Edge(pt2.index).coordOfOppositeVert(1);
				double b = sqrt(detaX * detaX + detaY * detaY);
				double leftAngle2 = acos((a * a + b * b - c * c) / (2.0 * a * b));
				double rightAngle2 = Edge(pt2.index).angleOpposite - leftAngle2;
				rightAngleSum += leftAngle2;
				leftAngleSum += rightAngle2;
			}
		}
		else
		{
			assert(!pt1.isVertex && !pt2.isVertex);
			int index1 = GetSubindexToVert(root, Edge(pt1.index).indexOfLeftVert);
			int index2 = GetSubindexToVert(root, Edge(pt1.index).indexOfRightVert);
			int index3 = GetSubindexToVert(root, Edge(pt2.index).indexOfLeftVert);
			int index4 = GetSubindexToVert(root, Edge(pt2.index).indexOfRightVert);
			if (index1 == index3)
			{
				double angleSum = 0;
				for (int i = 0; i < Neigh(root).size(); ++i)
					angleSum += Neigh(root)[i].second;
				if (abs(pt1.proportion - pt2.proportion) < 1e-2
					|| Edge(pt1.index).angleOpposite < 5.0 * M_PI / 180)
				{
					if (pt1.proportion > pt2.proportion)
					{
						leftAngleSum = (pt1.proportion - pt2.proportion) * Edge(pt1.index).angleOpposite;
						rightAngleSum = angleSum - leftAngleSum;
					}
					else
					{
						rightAngleSum = (pt2.proportion - pt1.proportion) * Edge(pt1.index).angleOpposite;
						leftAngleSum = angleSum - rightAngleSum;
					}
				}
				else if (pt1.proportion > pt2.proportion)
				{
					double c = (pt1.proportion - pt2.proportion) * Edge(pt1.index).length;
					double detaX1 = pt1.proportion * Edge(pt1.index).length - Edge(pt1.index).coordOfOppositeVert(0);
					double detaY1 = 0 - Edge(pt1.index).coordOfOppositeVert(1);
					double a = sqrt(detaX1 * detaX1 + detaY1 * detaY1);
					double detaX2 = pt2.proportion * Edge(pt2.index).length - Edge(pt2.index).coordOfOppositeVert(0);
					double detaY2 = 0 - Edge(pt2.index).coordOfOppositeVert(1);
					double b = sqrt(detaX2 * detaX2 + detaY2 * detaY2);
					leftAngleSum = acos((a * a + b * b - c * c) / (2 * a * b));
					rightAngleSum = angleSum - leftAngleSum;
				}
				else
				{
					double c = (pt1.proportion - pt2.proportion) * Edge(pt1.index).length;
					double detaX1 = pt1.proportion * Edge(pt1.index).length - Edge(pt1.index).coordOfOppositeVert(0);
					double detaY1 = 0 - Edge(pt1.index).coordOfOppositeVert(1);
					double a = sqrt(detaX1 * detaX1 + detaY1 * detaY1);
					double detaX2 = pt2.proportion * Edge(pt2.index).length - Edge(pt2.index).coordOfOppositeVert(0);
					double detaY2 = 0 - Edge(pt2.index).coordOfOppositeVert(1);
					double b = sqrt(detaX2 * detaX2 + detaY2 * detaY2);
					rightAngleSum = acos((a * a + b * b - c * c) / (2 * a * b));
					leftAngleSum = angleSum - rightAngleSum;
				}
			}
			else
			{
				double leftAngle1, rightAngle1, leftAngle2, rightAngle2;
				if (pt1.proportion < 1e-2
					|| pt1.proportion > 1 - 1e-2
					|| Edge(pt1.index).angleOpposite < 5.0 * M_PI / 180)
				{
					leftAngle1 = pt1.proportion * Edge(pt1.index).angleOpposite;
					rightAngle1 = (1 - pt1.proportion) * Edge(pt1.index).angleOpposite;
				}
				else
				{
					double c = Edge(pt1.index).length * pt1.proportion;
					double a = Edge(Edge(pt1.index).indexOfLeftEdge).length;
					double detaX = c - Edge(pt1.index).coordOfOppositeVert(0);
					double detaY = 0 - Edge(pt1.index).coordOfOppositeVert(1);
					double b = sqrt(detaX * detaX + detaY * detaY);
					leftAngle1 = acos((a * a + b * b - c * c) / (2.0 * a * b));
					rightAngle1 = Edge(pt1.index).angleOpposite - leftAngle1;
				}
				if (pt2.proportion < 1e-2
					|| pt2.proportion > 1 - 1e-2
					|| Edge(pt2.index).angleOpposite < 5.0 * M_PI / 180)
				{
					leftAngle2 = pt2.proportion * Edge(pt2.index).angleOpposite;
					rightAngle2 = (1 - pt2.proportion) * Edge(pt2.index).angleOpposite;
				}
				else
				{
					double c = Edge(pt2.index).length * pt2.proportion;
					double a = Edge(Edge(pt2.index).indexOfLeftEdge).length;
					double detaX = c - Edge(pt2.index).coordOfOppositeVert(0);
					double detaY = 0 - Edge(pt2.index).coordOfOppositeVert(1);
					double b = sqrt(detaX * detaX + detaY * detaY);
					leftAngle2 = acos((a * a + b * b - c * c) / (2.0 * a * b));
					rightAngle2 = Edge(pt2.index).angleOpposite - leftAngle2;
				}
				leftAngleSum = leftAngle1 + rightAngle2;
				rightAngleSum = rightAngle1 + leftAngle2;
				int index = index2;
				while (index != index3)
				{
					rightAngleSum += Neigh(root)[index].second;
					index = (index + 1) % Neigh(root).size();
				}
				index = index4;
				while (index != index1)
				{
					leftAngleSum += Neigh(root)[index].second;
					index = (index + 1) % Neigh(root).size();
				}
			}
		}
		return make_pair(leftAngleSum, rightAngleSum);
	}


	int CRichModel::GetNumOfValidDirectedEdges() const
	{
		return (int)m_Faces.size() * 3;
	}

	int CRichModel::GetNumOfTotalUndirectedEdges() const
	{
		return (int)m_Edges.size() / 2;
	}

	int CRichModel::GetNumOfGenera() const
	{
		return int(GetNumOfTotalUndirectedEdges() - (GetNumOfVerts() - m_nIsolatedVerts) - GetNumOfFaces() - GetNumOfBoundries()) / 2 + 1;
	}

	int CRichModel::GetNumOfComponents() const
	{
		return m_nComponents;
	}

	int CRichModel::GetNumOfBoundries() const
	{
		return m_nBoundries;
	}

	bool CRichModel::IsClosedModel() const
	{
		return GetNumOfValidDirectedEdges() == GetNumOfEdges();
	}

	int CRichModel::GetNumOfIsolated() const
	{
		return m_nIsolatedVerts;
	}

	int CRichModel::GetNumOfEdges() const
	{
		return (int)m_Edges.size();
	}

	bool CRichModel::isBoundaryVert(int index) const
	{
		return Neigh(index).empty() || IsStartEdge(Neigh(index).front().first);
	}

	bool CRichModel::IsStronglyConvexVert(int index) const
	{
		return m_FlagsForCheckingConvexVerts[index].first;
	}

	bool CRichModel::IsWeaklyConvexVert(int index) const
	{
		return m_FlagsForCheckingConvexVerts[index].second;
	}

	bool CRichModel::IsExtremeEdge(int edgeIndex) const
	{
		return Edge(edgeIndex).indexOfOppositeVert == -1;
	}

	bool CRichModel::IsStartEdge(int edgeIndex) const
	{
		return Edge(Edge(edgeIndex).indexOfReverseEdge).indexOfOppositeVert == -1;
	}

	const CRichModel::CEdge& CRichModel::Edge(int edgeIndex) const
	{
		return m_Edges[edgeIndex];
	}

	const vector<pair<int, double> >& CRichModel::Neigh(int root) const
	{
		return m_NeighsAndAngles[root];
	}

	double CRichModel::ProportionOnEdgeByImage(int edgeIndex, const Eigen::Vector2d& coord) const
	{
		double res = Edge(edgeIndex).coordOfOppositeVert(0) * coord(1) - Edge(edgeIndex).coordOfOppositeVert(1) * coord(0);
		return res / ((coord(1) - Edge(edgeIndex).coordOfOppositeVert(1)) * Edge(edgeIndex).length);
	}

	double CRichModel::ProportionOnEdgeByImage(int edgeIndex, double x1, double y1, double x2, double y2) const
	{
		double res = x1 * y2 - x2 * y1;
		return res / ((y2 - y1) * Edge(edgeIndex).length);
	}

	double CRichModel::ProportionOnLeftEdgeByImage(int edgeIndex, const Eigen::Vector2d &coord, double proportion) const
	{
		double xBalance = proportion * Edge(edgeIndex).length;
		double res = Edge(edgeIndex).coordOfOppositeVert(0) * coord(1) - Edge(edgeIndex).coordOfOppositeVert(1) * (coord(0) - xBalance);
		return xBalance * coord(1) / res;
	}

	double CRichModel::ProportionOnRightEdgeByImage(int edgeIndex, const Eigen::Vector2d &coord, double proportion) const
	{
		double part1 = Edge(edgeIndex).length * coord(1);
		double part2 = proportion * Edge(edgeIndex).length * Edge(edgeIndex).coordOfOppositeVert(1);
		double part3 = Edge(edgeIndex).coordOfOppositeVert(1) * coord(0) - Edge(edgeIndex).coordOfOppositeVert(0) * coord(1);
		return (part3 + proportion * part1 - part2) / (part3 + part1 - part2);
	}

	Eigen::Vector2d CRichModel::GetNew2DCoordinatesByRotatingAroundLeftChildEdge(int edgeIndex, const Eigen::Vector2d& input2DCoordinates) const
	{
		return Eigen::Vector2d(Edge(edgeIndex).matrixRotatedToLeftEdge(0) * input2DCoordinates(0) - Edge(edgeIndex).matrixRotatedToLeftEdge(1) * input2DCoordinates(1),
			Edge(edgeIndex).matrixRotatedToLeftEdge(1) * input2DCoordinates(0) + Edge(edgeIndex).matrixRotatedToLeftEdge(0) * input2DCoordinates(1));
	}

	Eigen::Vector2d CRichModel::GetNew2DCoordinatesByRotatingAroundRightChildEdge(int edgeIndex, const Eigen::Vector2d& input2DCoordinates) const
	{
		int reverseEdge = Edge(Edge(edgeIndex).indexOfRightEdge).indexOfReverseEdge;
		Eigen::Vector2d coordOfLeftEnd = GetNew2DCoordinatesByReversingCurrentEdge(reverseEdge, Edge(reverseEdge).coordOfOppositeVert);
		return Eigen::Vector2d(Edge(edgeIndex).matrixRotatedToRightEdge(0) * input2DCoordinates(0) - Edge(edgeIndex).matrixRotatedToRightEdge(1) * input2DCoordinates(1) + coordOfLeftEnd(0),
			Edge(edgeIndex).matrixRotatedToRightEdge(1) * input2DCoordinates(0) + Edge(edgeIndex).matrixRotatedToRightEdge(0) * input2DCoordinates(1) + coordOfLeftEnd(1));
	}

	Eigen::Vector2d CRichModel::GetNew2DCoordinatesByReversingCurrentEdge(int edgeIndex, const Eigen::Vector2d& input2DCoordinates) const
	{
		return Eigen::Vector2d(Edge(edgeIndex).length - input2DCoordinates(0), -input2DCoordinates(1));
	}

	int CRichModel::GetSubindexToVert(int root, int neigh) const
	{
		for (int i = 0; i < (int)Neigh(root).size(); ++i)
		{
			if (Edge(Neigh(root)[i].first).indexOfRightVert == neigh)
				return i;
		}
		return -1;
	}

	double CRichModel::DistanceToOppositeAngle(int edgeIndex, const Eigen::Vector2d& coord) const
	{
		double detaX = coord(0) - Edge(edgeIndex).coordOfOppositeVert(0);
		double detaY = coord(1) - Edge(edgeIndex).coordOfOppositeVert(1);
		return sqrt(detaX * detaX + detaY * detaY);
	}

	double CRichModel::DistanceToLeftVert(int edgeIndex, const Eigen::Vector2d& coord) const
	{
		double detaX = coord(0);
		double detaY = coord(1);
		return sqrt(detaX * detaX + detaY * detaY);
	}

	double CRichModel::DistanceToRightVert(int edgeIndex, const Eigen::Vector2d& coord) const
	{
		double detaX = coord(0) - Edge(edgeIndex).length;
		double detaY = coord(1);
		return sqrt(detaX * detaX + detaY * detaY);
	}

	int CRichModel::GetEdgeIndexFromTwoVertices(int leftVert, int rightVert) const
	{
		int subIndex = GetSubindexToVert(leftVert, rightVert);
		assert(subIndex != -1);
		return Neigh(leftVert)[subIndex].first;
	}

	int CRichModel::SplitEdge(const EdgePoint& ep)
	{
		CPoint3D newVert = (1 - ep.proportion) * m_Verts[Edge(ep.index).indexOfLeftVert]
			+ ep.proportion * m_Verts[Edge(ep.index).indexOfRightVert];
		m_Verts.push_back(newVert);
		int vertID = m_Verts.size() - 1;
		int edgeIndex = ep.index;
		m_UselessFaces.insert(Edge(edgeIndex).indexOfFrontFace);
		int reverseFaceIndex = Edge(Edge(edgeIndex).indexOfReverseEdge).indexOfFrontFace;
		if (reverseFaceIndex != -1)
			m_UselessFaces.insert(reverseFaceIndex);

		m_UselessEdges.insert(edgeIndex);
		if (reverseFaceIndex != -1)
			m_UselessEdges.insert(Edge(edgeIndex).indexOfReverseEdge);
		m_Faces.push_back(CBaseModel::CFace(Edge(edgeIndex).indexOfLeftVert, vertID, Edge(edgeIndex).indexOfOppositeVert));
		int f1 = (int)m_Faces.size() - 1;
		m_Faces.push_back(CBaseModel::CFace(vertID, m_Edges[edgeIndex].indexOfRightVert, m_Edges[edgeIndex].indexOfOppositeVert));
		int f2 = (int)m_Faces.size() - 1;
		if (reverseFaceIndex != -1)
			m_Faces.push_back(CBaseModel::CFace(m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfLeftVert, vertID, m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfOppositeVert));
		int f3 = (int)m_Faces.size() - 1;
		if (reverseFaceIndex != -1)
			m_Faces.push_back(CBaseModel::CFace(vertID, m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfRightVert, m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfOppositeVert));
		int f4 = (int)m_Faces.size() - 1;
		if (reverseFaceIndex == -1)
			f3 = f4 = -1;
		m_Edges.push_back(CRichModel::CEdge());
		int e1 = (int)m_Edges.size() - 1;
		m_Edges.push_back(CRichModel::CEdge());
		int e2 = (int)m_Edges.size() - 1;
		m_Edges.push_back(CRichModel::CEdge());
		int e3 = (int)m_Edges.size() - 1;
		m_Edges.push_back(CRichModel::CEdge());
		int e4 = (int)m_Edges.size() - 1;
		m_Edges.push_back(CRichModel::CEdge());
		int e5 = (int)m_Edges.size() - 1;
		if (reverseFaceIndex != -1)
			m_Edges.push_back(CRichModel::CEdge());
		int e6 = (int)m_Edges.size() - 1;
		if (reverseFaceIndex != -1)
			m_Edges.push_back(CRichModel::CEdge());
		int e7 = (int)m_Edges.size() - 1;
		m_Edges.push_back(CRichModel::CEdge());
		int e8 = (int)m_Edges.size() - 1;
		if (reverseFaceIndex == -1)
			e6 = e7 = -1;

		m_Edges[e1].indexOfFrontFace = f1;
		m_Edges[e1].indexOfLeftEdge = m_Edges[edgeIndex].indexOfLeftEdge;
		m_Edges[e1].indexOfLeftVert = m_Edges[edgeIndex].indexOfLeftVert;
		m_Edges[e1].indexOfOppositeVert = m_Edges[edgeIndex].indexOfOppositeVert;
		m_Edges[e1].indexOfReverseEdge = e8;
		m_Edges[e1].indexOfRightEdge = e3;
		m_Edges[e1].indexOfRightVert = vertID;

		m_Edges[e4].indexOfFrontFace = f2;
		m_Edges[e4].indexOfLeftEdge = e2;
		m_Edges[e4].indexOfLeftVert = vertID;
		m_Edges[e4].indexOfOppositeVert = m_Edges[edgeIndex].indexOfOppositeVert;
		m_Edges[e4].indexOfReverseEdge = e5;
		m_Edges[e4].indexOfRightEdge = m_Edges[edgeIndex].indexOfRightEdge;
		m_Edges[e4].indexOfRightVert = m_Edges[edgeIndex].indexOfRightVert;

		m_Edges[e5].indexOfFrontFace = f3;
		m_Edges[e5].indexOfLeftEdge = m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfLeftEdge;
		m_Edges[e5].indexOfLeftVert = m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfLeftVert;
		m_Edges[e5].indexOfOppositeVert = m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfOppositeVert;
		m_Edges[e5].indexOfReverseEdge = e4;
		m_Edges[e5].indexOfRightEdge = e7;
		m_Edges[e5].indexOfRightVert = vertID;

		m_Edges[e8].indexOfFrontFace = f4;
		m_Edges[e8].indexOfLeftEdge = e6;
		m_Edges[e8].indexOfLeftVert = vertID;
		m_Edges[e8].indexOfOppositeVert = m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfOppositeVert;
		m_Edges[e8].indexOfReverseEdge = e1;
		m_Edges[e8].indexOfRightEdge = m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfRightEdge;
		m_Edges[e8].indexOfRightVert = m_Edges[m_Edges[edgeIndex].indexOfReverseEdge].indexOfRightVert;

		m_Edges[e2].indexOfFrontFace = f1;
		m_Edges[e2].indexOfLeftEdge = m_Edges[e1].indexOfReverseEdge;
		m_Edges[e2].indexOfLeftVert = m_Edges[e1].indexOfRightVert;
		m_Edges[e2].indexOfOppositeVert = m_Edges[e1].indexOfLeftVert;
		m_Edges[e2].indexOfReverseEdge = e3;
		m_Edges[e2].indexOfRightEdge = m_Edges[e1].indexOfLeftEdge;
		m_Edges[e2].indexOfRightVert = m_Edges[e1].indexOfOppositeVert;

		m_Edges[e3].indexOfFrontFace = f2;
		m_Edges[e3].indexOfLeftEdge = m_Edges[e4].indexOfRightEdge;
		m_Edges[e3].indexOfLeftVert = m_Edges[e4].indexOfOppositeVert;
		m_Edges[e3].indexOfOppositeVert = m_Edges[e4].indexOfRightVert;
		m_Edges[e3].indexOfReverseEdge = e2;
		m_Edges[e3].indexOfRightEdge = m_Edges[e4].indexOfReverseEdge;
		m_Edges[e3].indexOfRightVert = m_Edges[e4].indexOfLeftVert;

		if (reverseFaceIndex != -1)
		{
			m_Edges[e6].indexOfFrontFace = f3;
			m_Edges[e6].indexOfLeftEdge = m_Edges[e5].indexOfReverseEdge;
			m_Edges[e6].indexOfLeftVert = m_Edges[e5].indexOfRightVert;
			m_Edges[e6].indexOfOppositeVert = m_Edges[e5].indexOfLeftVert;
			m_Edges[e6].indexOfReverseEdge = e7;
			m_Edges[e6].indexOfRightEdge = m_Edges[e5].indexOfLeftEdge;
			m_Edges[e6].indexOfRightVert = m_Edges[e5].indexOfOppositeVert;

			m_Edges[e7].indexOfFrontFace = f4;
			m_Edges[e7].indexOfLeftEdge = m_Edges[e8].indexOfRightEdge;
			m_Edges[e7].indexOfLeftVert = m_Edges[e8].indexOfOppositeVert;
			m_Edges[e7].indexOfOppositeVert = m_Edges[e8].indexOfRightVert;
			m_Edges[e7].indexOfReverseEdge = e6;
			m_Edges[e7].indexOfRightEdge = m_Edges[e8].indexOfReverseEdge;
			m_Edges[e7].indexOfRightVert = m_Edges[e8].indexOfLeftVert;
		}
		int newEdgeCnt = 8;
		if (reverseFaceIndex == -1)
			newEdgeCnt = 6;
		//for (int i = 0; i < newEdgeCnt; ++i)
		//{
		//	extraEdgePool[make_pair(newEdges[newEdges.size() - 1 - i].indexOfLeftVert, newEdges[newEdges.size() - 1 - i].indexOfRightVert)] = newEdges.size() - 1 - i;
		//}

		int leftEdge = m_Edges[m_Edges[e1].indexOfLeftEdge].indexOfReverseEdge;
		m_Edges[leftEdge].indexOfFrontFace = f1;
		m_Edges[leftEdge].indexOfLeftEdge = e3;
		m_Edges[leftEdge].indexOfLeftVert;
		m_Edges[leftEdge].indexOfOppositeVert = m_Edges[e1].indexOfRightVert;
		m_Edges[leftEdge].indexOfReverseEdge;
		m_Edges[leftEdge].indexOfRightEdge = m_Edges[e1].indexOfReverseEdge;
		m_Edges[leftEdge].indexOfRightVert;

		int rightEdge = m_Edges[m_Edges[e4].indexOfRightEdge].indexOfReverseEdge;
		m_Edges[rightEdge].indexOfFrontFace = f2;
		m_Edges[rightEdge].indexOfLeftEdge = m_Edges[e4].indexOfReverseEdge;
		m_Edges[rightEdge].indexOfLeftVert;
		m_Edges[rightEdge].indexOfOppositeVert = m_Edges[e4].indexOfLeftVert;
		m_Edges[rightEdge].indexOfReverseEdge;
		m_Edges[rightEdge].indexOfRightEdge = e2;
		m_Edges[rightEdge].indexOfRightVert;

		if (reverseFaceIndex != -1)
		{
			leftEdge = m_Edges[m_Edges[e5].indexOfLeftEdge].indexOfReverseEdge;
			m_Edges[leftEdge].indexOfFrontFace = f3;
			m_Edges[leftEdge].indexOfLeftEdge = e7;
			m_Edges[leftEdge].indexOfLeftVert;
			m_Edges[leftEdge].indexOfOppositeVert = m_Edges[e5].indexOfRightVert;
			m_Edges[leftEdge].indexOfReverseEdge;
			m_Edges[leftEdge].indexOfRightEdge = m_Edges[e5].indexOfReverseEdge;
			m_Edges[leftEdge].indexOfRightVert;

			rightEdge = m_Edges[m_Edges[e8].indexOfRightEdge].indexOfReverseEdge;
			m_Edges[rightEdge].indexOfFrontFace = f4;
			m_Edges[rightEdge].indexOfLeftEdge = m_Edges[e8].indexOfReverseEdge;
			m_Edges[rightEdge].indexOfLeftVert;
			m_Edges[rightEdge].indexOfOppositeVert = m_Edges[e8].indexOfLeftVert;
			m_Edges[rightEdge].indexOfReverseEdge;
			m_Edges[rightEdge].indexOfRightEdge = e6;
			m_Edges[rightEdge].indexOfRightVert;
		}
		return vertID;
	}

	void CRichModel::SavePathToObj(const vector<EdgePoint>& pl, const string& filename) const
	{
		ofstream out(filename.c_str());
		//filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1);
		//filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1)
		out << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
		if (!pl.empty())
		{
			for (int i = 0; i < pl.size(); ++i)
			{
				CPoint3D pt = pl[i].GetShiftPoint(*this);
				out << "v " << pt.x << " " << pt.y << " " << pt.z << endl;
			}

			out << "l ";
			for (int i = 0; i < pl.size(); ++i)
			{
				out << i + 1 << " ";
			}
			out << endl;
		}
		out.close();
	}


	//void CRichModel::SavePathToObj(const vector<CPoint3D>& pl, const string& filename) const
	//{
	//	ofstream out(filename.c_str());
	//	out << "g 3D_Curve" << endl;
	//	if (!pl.empty())
	//	{	
	//		for (int i = 0; i < pl.size(); ++i)
	//		{
	//			CPoint3D pt = pl[i];
	//			out << "v " << pt.x << " " << pt.y << " " << pt.z << endl;
	//		}
	//
	//		out << "l ";
	//		for (int i = 0; i < pl.size(); ++i)
	//		{
	//			out << i + 1 << " ";
	//		}
	//		out << endl;
	//	}
	//	out.close();
	//}

	

	void CRichModel::SaveIsolineToObj(const vector<EdgePoint>& isoline, const string& filename) const
	{
		ofstream out(filename.c_str());
		out << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
		if (!isoline.empty())
		{
			for (int i = 0; i < isoline.size(); ++i)
			{
				CPoint3D pt = isoline[i].GetShiftPoint(*this);
				out << "v " << pt.x << " " << pt.y << " " << pt.z << endl;
			}
			out << "l ";
			for (int i = 0; i < isoline.size(); ++i)
			{
				out << i + 1 << " ";
			}
			out << "1\n";
		}
		out.close();
	}

	void CRichModel::SaveIsolineToObj(const vector<vector<EdgePoint>>& isoline, const string& filename) const
	{
		ofstream out(filename.c_str());
		out << "g " << filename.substr(filename.rfind("\\") + 1, filename.rfind('.') - filename.rfind("\\") - 1) << endl;
		int base = 0;
		for (int nth = 0; nth < isoline.size(); ++nth)
		{
			for (int i = 0; i < isoline[nth].size(); ++i)
			{
				CPoint3D pt = isoline[nth][i].GetShiftPoint(*this);
				out << "v " << pt.x << " " << pt.y << " " << pt.z << endl;
			}
			out << "l ";
			for (int i = 0; i < isoline[nth].size(); ++i)
			{
				out << i + 1 + base << " ";
			}
			out << base + 1 << "\n";
			base += isoline[nth].size();
		}
		out.close();
	}
	pair<CBaseModel, CBaseModel> CRichModel::SplitBasedOnScalarField_into_Two_Models(const vector<double>& scalarField,
		double val) const
	{
		CRichModel clone(*this);
		vector<EdgePoint> eps;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (Edge(i).indexOfLeftVert >
				Edge(i).indexOfRightVert)
				continue;
			if (scalarField[Edge(i).indexOfLeftVert] >= val
				== scalarField[Edge(i).indexOfRightVert] >= val)
				continue;
			double prop = (val - scalarField[Edge(i).indexOfLeftVert])
				/ (scalarField[Edge(i).indexOfRightVert] - scalarField[Edge(i).indexOfLeftVert]);
			if (prop < 1e-4 || prop > 1 - 1e-4)
				continue;
			eps.push_back(EdgePoint(i, prop));
		}
		int oldVertNum = GetNumOfVerts();
		for (int i = 0; i < eps.size(); ++i)
			clone.SplitEdge(eps[i]);
		vector<CPoint3D> vertList(clone.m_Verts);
		vector<CFace> faceList_small;
		vector<CFace> faceList_large;
		for (int i = 0; i < clone.GetNumOfFaces(); ++i)
		{
			if (clone.m_UselessFaces.find(i) != clone.m_UselessFaces.end())
				continue;
			int v1 = clone.Face(i)[0];
			int v2 = clone.Face(i)[1];
			int v3 = clone.Face(i)[2];
			int cnt(0);
			double average(0);
			if (v1 < oldVertNum)
			{
				average += scalarField[v1];
				cnt++;
			}
			if (v2 < oldVertNum)
			{
				average += scalarField[v2];
				cnt++;
			}
			if (v3 < oldVertNum)
			{
				average += scalarField[v3];
				cnt++;
			}
			average /= cnt;
			if (average < val)
			{
				faceList_small.push_back(clone.Face(i));
			}
			else
			{
				faceList_large.push_back(clone.Face(i));
			}
		}
		CBaseModel model_small(vertList, faceList_small);
		model_small.RemoveUnreferencedVertices();
		CBaseModel model_large(vertList, faceList_large);
		model_large.RemoveUnreferencedVertices();
		return make_pair(model_small, model_large);
	}
	void CRichModel::SplitBasedOnScalarField_Update(const vector<double>& scalarField,
		double val)
	{
		vector<EdgePoint> eps;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (Edge(i).indexOfLeftVert >
				Edge(i).indexOfRightVert)
				continue;
			if (scalarField[Edge(i).indexOfLeftVert] >= val
				== scalarField[Edge(i).indexOfRightVert] >= val)
				continue;
			double prop = (val - scalarField[Edge(i).indexOfLeftVert])
				/ (scalarField[Edge(i).indexOfRightVert] - scalarField[Edge(i).indexOfLeftVert]);
			if (prop < 1e-4 || prop > 1 - 1e-4)
				continue;
			eps.push_back(EdgePoint(i, prop));
		}
		for (int i = 0; i < eps.size(); ++i)
			SplitEdge(eps[i]);
		vector<CFace> faceList;
		for (int i = 0; i < GetNumOfFaces(); ++i)
		{
			if (m_UselessFaces.find(i) != m_UselessFaces.end())
				continue;
			faceList.push_back(Face(i));			
		}
		swap(m_Faces, faceList);
		PreprocessVertsAndFacesIntoBaseModel();
		PreprocessBaseModelIntoRichModel();
	}

	void CRichModel::SplitBasedOnScalarField_into_Two_Models(const vector<double>& scalarField,
		double val,
		const string& fileWithLargerScalars,
		const string& fileWithSmallerScalars) const
	{
		CRichModel clone(*this);
		vector<EdgePoint> eps;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (Edge(i).indexOfLeftVert >
				Edge(i).indexOfRightVert)
				continue;
			if (scalarField[Edge(i).indexOfLeftVert] >= val
				== scalarField[Edge(i).indexOfRightVert] >= val)
				continue;
			double prop = (val - scalarField[Edge(i).indexOfLeftVert])
				/ (scalarField[Edge(i).indexOfRightVert] - scalarField[Edge(i).indexOfLeftVert]);
			if (prop < 1e-4 || prop > 1 - 1e-4)
				continue;
			eps.push_back(EdgePoint(i, prop));
		}
		int oldVertNum = GetNumOfVerts();
		for (int i = 0; i < eps.size(); ++i)
			clone.SplitEdge(eps[i]);

		ofstream outLarge(fileWithLargerScalars.c_str());
		ofstream outSmall(fileWithSmallerScalars.c_str());
		for (int i = 0; i < clone.GetNumOfVerts(); ++i)
		{
			outLarge << "v " << clone.Vert(i).x
				<< " " << clone.Vert(i).y
				<< " " << clone.Vert(i).z << endl;
			outSmall << "v " << clone.Vert(i).x
				<< " " << clone.Vert(i).y
				<< " " << clone.Vert(i).z << endl;
		}

		for (int i = 0; i < clone.GetNumOfFaces(); ++i)
		{
			if (clone.m_UselessFaces.find(i) != clone.m_UselessFaces.end())
				continue;
			int v1 = clone.Face(i)[0];
			int v2 = clone.Face(i)[1];
			int v3 = clone.Face(i)[2];
			int cnt(0);
			double average(0);
			if (v1 < oldVertNum)
			{
				average += scalarField[v1];
				cnt++;
			}
			if (v2 < oldVertNum)
			{
				average += scalarField[v2];
				cnt++;
			}
			if (v3 < oldVertNum)
			{
				average += scalarField[v3];
				cnt++;
			}
			average /= cnt;
			if (average < val)
			{
				outSmall << "f " << clone.Face(i)[0] + 1 << " " << clone.Face(i)[1] + 1 << " " << clone.Face(i)[2] + 1 << endl;
			}
			else
			{
				outLarge << "f " << clone.Face(i)[0] + 1 << " " << clone.Face(i)[1] + 1 << " " << clone.Face(i)[2] + 1 << endl;
			}
		}

		outLarge.close();
		outSmall.close();
	}

	int CRichModel::IntersectQuery(int faceID, const pair<EdgePoint, EdgePoint>& seg1, const pair<EdgePoint, EdgePoint>& seg2, EdgePoint& intersection) const
	{
		set<pair<double, string>> sortingSet;
		sortingSet.insert(make_pair(seg1.first.GetNumbering(*this, faceID), "seg1"));
		sortingSet.insert(make_pair(seg1.second.GetNumbering(*this, faceID), "seg1"));
		sortingSet.insert(make_pair(seg2.first.GetNumbering(*this, faceID), "seg2"));
		sortingSet.insert(make_pair(seg2.second.GetNumbering(*this, faceID), "seg2"));
		vector<pair<double, string>> sortingVec(sortingSet.begin(), sortingSet.end());
		for (int i = 0; i < 4; ++i)
		{
			int nxt = (i + 1) % 4;
			if (sortingVec[i].first == sortingVec[nxt].first
				&& sortingVec[i].second != sortingVec[nxt].second)
			{
				if (seg1.first == seg2.first || seg1.first == seg2.second)
					intersection = seg1.first;
				else
					intersection = seg1.second;
				return -1;
			}
		}
		for (int i = 0; i < 4; ++i)
		{
			int nxt = (i + 1) % 4;
			if (sortingVec[i].second == sortingVec[nxt].second)
				return 0;
		}
		return 1;
	}

	double CRichModel::GetMaxEdgeLength() const
	{
		return m_maxEdgeLength;
	}

	void CRichModel::SplitEdgeSet(const set<int> & edgeSet, double tolerance)
	{
		auto edgeSet_clr = edgeSet;
		edgeSet_clr.clear();
		for (auto e : edgeSet)
		{
			auto e_reverse = Edge(e).indexOfReverseEdge;
			if (edgeSet.find(e_reverse) == edgeSet.end() || Edge(e).indexOfLeftVert < Edge(e).indexOfRightVert)
				edgeSet_clr.insert(e);
		}
		for (auto e : edgeSet_clr)
		{
			auto e_reverse = Edge(e).indexOfReverseEdge;
			if (Edge(e).indexOfFrontFace != -1)
			{
				m_UselessFaces.insert(Edge(e).indexOfFrontFace);
			}
			if (Edge(e_reverse).indexOfFrontFace != -1)
			{
				m_UselessFaces.insert(Edge(e_reverse).indexOfFrontFace);
			}
			int oldVertSize = m_Verts.size();
			int segNum = int(Edge(e).length / tolerance + 1);
			double gap = Edge(e).length / segNum;
			for (int j = 1; j < segNum; ++j)
			{
				m_Verts.push_back(Vert(Edge(e).indexOfLeftVert) + j / (double)segNum * (Vert(Edge(e).indexOfRightVert) - Vert(Edge(e).indexOfLeftVert)));
			}
			for (int j = 1; j <= segNum; ++j)
			{
				int vID1 = j + oldVertSize - 2;
				int vID2 = j + oldVertSize - 1;
				if (j == 1)
					vID1 = Edge(e).indexOfLeftVert;
				if (j == segNum)
					vID2 = Edge(e).indexOfRightVert;
				if (Edge(e).indexOfOppositeVert != -1)
					m_Faces.push_back(CFace(vID1, vID2, Edge(e).indexOfOppositeVert));
				if (Edge(e_reverse).indexOfOppositeVert != -1)
					m_Faces.push_back(CFace(vID2, vID1, Edge(e_reverse).indexOfOppositeVert));
			}
		}
	}

	vector<vector<EdgePoint>> CRichModel::GetIsoLine(const vector<double>& scalarField,
		double val) const
	{
		vector<vector<EdgePoint>> result;
		map<int, EdgePoint> intersections;
		set<int> vertices;
		for (int i = 0; i < GetNumOfEdges(); ++i)
		{
			if (Edge(i).indexOfLeftVert >
				Edge(i).indexOfRightVert)
				continue;

			if (scalarField[Edge(i).indexOfLeftVert] >= val
				== scalarField[Edge(i).indexOfRightVert] >= val)
				continue;
			double prop = (val - scalarField[Edge(i).indexOfLeftVert])
				/ (scalarField[Edge(i).indexOfRightVert] - scalarField[Edge(i).indexOfLeftVert]);
			if (prop < 2e-2)
			{
				vertices.insert(Edge(i).indexOfLeftVert);
			}
			else if (prop > 1 - 2e-2)
			{
				vertices.insert(Edge(i).indexOfRightVert);
			}
			else
			{
				intersections[i] = EdgePoint(i, prop);
				intersections[Edge(i).indexOfReverseEdge] = EdgePoint(Edge(i).indexOfReverseEdge, 1 - prop);
			}
		}
		set<int> badEdgeIndices;
		for (auto mypair : intersections)
		{
			if (vertices.find(Edge(mypair.second.index).indexOfLeftVert) != vertices.end()
				||
				vertices.find(Edge(mypair.second.index).indexOfRightVert) != vertices.end())
			{
				badEdgeIndices.insert(mypair.second.index);
			}
		}
		for (auto index : badEdgeIndices)
			intersections.erase(index);
		while (!intersections.empty())
		{
			EdgePoint top = intersections.begin()->second;

			if (scalarField[Edge(top.index).indexOfLeftVert] < scalarField[Edge(top.index).indexOfRightVert])
				top = intersections[Edge(top.index).indexOfReverseEdge];
			intersections.erase(top.index);
			//cerr << "leftv = " << Edge(top.index).indexOfLeftVert << endl;
			//cerr << "rightv = " << Edge(top.index).indexOfRightVert << endl;
			intersections.erase(Edge(top.index).indexOfReverseEdge);
			vector<EdgePoint> loop;
			EdgePoint ep = top;
			bool flag = true;

			do
			{
				loop.push_back(ep);
				if (!ep.isVertex)
				{
					if (Edge(ep.index).indexOfLeftEdge == top.index
						|| Edge(ep.index).indexOfRightEdge == top.index)
						break;
					if (vertices.find(Edge(ep.index).indexOfOppositeVert) != vertices.end())
					{
						//cerr << "vert " << Edge(ep.index).indexOfOppositeVert << " to be deleted\n";
						vertices.erase(Edge(ep.index).indexOfOppositeVert);
						ep = EdgePoint(Edge(ep.index).indexOfOppositeVert);
					}
					else if (intersections.find(Edge(ep.index).indexOfLeftEdge) != intersections.end())
					{
						ep = intersections[Edge(ep.index).indexOfLeftEdge];
						intersections.erase(ep.index);
						intersections.erase(Edge(ep.index).indexOfReverseEdge);

					}
					else if (intersections.find(Edge(ep.index).indexOfRightEdge) != intersections.end())
					{
						ep = intersections[Edge(ep.index).indexOfRightEdge];
						intersections.erase(ep.index);
						intersections.erase(Edge(ep.index).indexOfReverseEdge);
					}
					else
					{
						ep = EdgePoint(Edge(ep.index).indexOfOppositeVert);
						//flag = false;
						//cerr << "face ahead: " << Edge(ep.index).indexOfFrontFace << endl;
						//cerr << "face past: " << Edge(Edge(ep.index).indexOfReverseEdge).indexOfFrontFace << endl;
						//cerr << "left v: " << Edge(ep.index).indexOfLeftVert << endl;
						//cerr << "right v: " << Edge(ep.index).indexOfRightVert << endl;
						//cerr << "oppo v: " << Edge(ep.index).indexOfOppositeVert << endl;
						//SavePathToObj(loop, "loop1.obj");
						//break;
					}
				}
				else
				{
					bool fFound = false;
					for (int k = 0; k < Neigh(ep.index).size(); ++k)
					{
						int neigh = Edge(Neigh(ep.index)[k].first).indexOfRightVert;
						if (vertices.find(neigh) != vertices.end())
						{
							fFound = true;
							vertices.erase(neigh);
							ep = EdgePoint(neigh);
							break;
						}
					}
					if (!fFound)
					{
						fFound = false;
						bool done = false;
						for (int k = 0; k < Neigh(ep.index).size(); ++k)
						{
							int nextEdge = Edge(Neigh(ep.index)[k].first).indexOfRightEdge;
							//cerr << "leftV: " << Edge(nextEdge).indexOfLeftVert << endl;
							//cerr << "rightV: " << Edge(nextEdge).indexOfRightVert << endl;
							//cerr << "--leftV: " << Edge(top.index).indexOfLeftVert << endl;
							//cerr << "--rightV: " << Edge(top.index).indexOfRightVert << endl;
							if (nextEdge == top.index)
							{
								//should skip out
								done = true;
								fFound = true;
								break;
							}
							if (intersections.find(nextEdge) != intersections.end())
							{
								ep = intersections[nextEdge];
								intersections.erase(nextEdge);
								intersections.erase(Edge(nextEdge).indexOfReverseEdge);
								fFound = true;
								break;
							}
						}
						if (done)
							break;
					}
					if (!fFound)
					{
						flag = false;
						//SavePathToObj(loop, "loop2.obj");
						break;
					}
				}
			} while (true);
			if (flag)
			{
				map<int, int> vIDs;
				int v2 = -1;
				int pos2 = -1;
				for (int kk = 0; kk < loop.size(); ++kk)
				{
					auto ep = loop[kk];
					if (ep.isVertex)
					{
						if (vIDs.find(ep.index) == vIDs.end())
							vIDs[ep.index] = 1;
						else
						{
							vIDs[ep.index] = vIDs[ep.index] + 1;
							v2 = ep.index;
							pos2 = kk;
						}
					}
				}
				if (v2 != -1)
				{
					bool fFirst = true;
					vector<EdgePoint> track1;
					vector<EdgePoint> track2;
					vector<EdgePoint>* ptrack = &track1;
					for (int kk = 0; kk < loop.size(); ++kk)
					{
						auto ep = loop[kk];
						ptrack->push_back(ep);
						if (ep.isVertex && ep == v2 && fFirst)
						{
							fFirst = false;
							ptrack = &track2;
						}
						else if (ep.isVertex && ep == v2 && !fFirst)
						{
							ptrack = &track1;
						}
					}
					if (track1.size() > track2.size())
						result.push_back(track1);
					else
						result.push_back(track2);
				}
				else
					result.push_back(loop);
			}
		}
		return result;
	}

	void CRichModel::ExtractBoundary()
	{
		m_BoundaryCurves.clear();
		m_OnWhichBoundary.clear();
		set<int> boundaryVerts;
		for (int i = 0; i < GetNumOfVerts(); ++i)
		{
			if (isBoundaryVert(i))
				boundaryVerts.insert(i);
		}

		while (!boundaryVerts.empty())
		{
			vector<int> cluster;
			int head = *boundaryVerts.begin();
			cluster.push_back(head);
			boundaryVerts.erase(boundaryVerts.begin());
			if (Neigh(cluster.back()).empty())
				continue;
			int candidate = Edge(Neigh(cluster.back()).back().first).indexOfRightVert;
			while (candidate != cluster.front())
			{
				cluster.push_back(candidate);
				boundaryVerts.erase(candidate);
				candidate = Edge(Neigh(cluster.back()).back().first).indexOfRightVert;
			}
			m_BoundaryCurves.push_back(cluster);
			for (auto v : cluster)
			{
				m_OnWhichBoundary[v] = m_BoundaryCurves.size() - 1;
			}
		}
	}

	vector<double> CRichModel::ReadScalarField(const char* file) const
	{
		ifstream in(file);
		double dis;
		vector<double> field;
		while (in >> dis)
			field.push_back(dis);
		in.close();
		return field;
	}

	Eigen::Vector2d CRichModel::Get2DCoord(int face, CPoint3D pt) const
	{
		CPoint3D origin = Vert(Face(face)[0]);
		CPoint3D dir1 = Vert(Face(face)[1]) - Vert(Face(face)[0]);
		dir1.Normalize();
		CPoint3D faceNormal = VectorCross(Vert(Face(face)[0]), Vert(Face(face)[1]), Vert(Face(face)[2]));
		faceNormal.Normalize();
		CPoint3D dir2 = faceNormal * dir1;
		dir2.Normalize();
		return Eigen::Vector2d(((pt - origin) ^ dir1), ((pt - origin) ^ dir2));
	}

	CPoint3D CRichModel::Get3DCoord(int face, Eigen::Vector2d pt) const
	{
		CPoint3D origin = Vert(Face(face)[0]);
		CPoint3D dir1 = Vert(Face(face)[1]) - Vert(Face(face)[0]);
		dir1.Normalize();
		CPoint3D faceNormal = VectorCross(Vert(Face(face)[0]), Vert(Face(face)[1]), Vert(Face(face)[2]));
		faceNormal.Normalize();
		CPoint3D dir2 = faceNormal * dir1;
		dir2.Normalize();
		return origin + pt(0) * dir1 + pt(1) * dir2;
	}
	void CRichModel::RemoveUnreferencedVertices()
	{
		CBaseModel::RemoveUnreferencedVertices();
		PreprocessBaseModelIntoRichModel();
	}
}