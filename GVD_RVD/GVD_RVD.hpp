#pragma once
#include <Eigen\dense>
#include <iostream>
#include <sstream>
#include <fstream>
#include "tsl/robin_set.h"
#include "tsl/robin_map.h"
#include <vector>
#include <set>
#include <omp.h>
#include <queue>
#include<math.h>
#include <map>
#include <algorithm>
#include "windows.h"
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
//#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Aff_transformation_3.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel inexact_Kernel;
typedef CGAL::Polyhedron_3<inexact_Kernel>  Polyhedron;
typedef CGAL::Polygon_2<inexact_Kernel> Polygon_2;
typedef CGAL::Nef_polyhedron_3<inexact_Kernel>  Nef_polyhedron;
typedef inexact_Kernel::Point_3 Point_3;
typedef inexact_Kernel::Point_2 Point_2;
typedef inexact_Kernel::Vector_3 Vector_3;
typedef inexact_Kernel::Plane_3 Plane_3;
typedef CGAL::Aff_transformation_3<inexact_Kernel> Aff_transformation_3;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
using namespace std;

#include "..\\Geodesic\Xin_Wang.h"
#if defined(_DEBUG) && defined(_WIN64)
#pragma comment(lib, "..\\lib\\Debug\\Geodesic.lib")
#endif
#if !defined(_DEBUG) && defined(_WIN64)
#pragma comment(lib, "..\\lib\\Release\\Geodesic.lib")
#endif

#include "..\\Model3D\RichModel.h"
#if defined(_DEBUG) && defined(_WIN64)
#pragma comment(lib, "..\\Debug\\lib\\Model3D.lib")
#endif
#if !defined(_DEBUG) && defined(_WIN64)
#pragma comment(lib, "..\\Release\\lib\\Model3D.lib")
#endif

#include <Eigen/dense>
#include <unordered_set>

namespace GVD
{
	struct OverPropagation : public Geodesic::CXin_Wang
	{
		//map<int, set<int>> m_destinationVerticesOfAncestors;
		map<int, set<int>> m_candidateCoveringFacesForEachAncestor;
		OverPropagation(const Model3D::CRichModel& model, const map<int, double>& sources)
			: Geodesic::CXin_Wang(model, sources) {
			for (auto mypair : sources)
			{
				for (int k = 0; k < model.Neigh(mypair.first).size(); ++k)
				{
					m_candidateCoveringFacesForEachAncestor[mypair.first].insert(model.Edge(model.Neigh(mypair.first)[k].first).indexOfFrontFace);
				}
			}
		}
		OverPropagation(const Model3D::CRichModel& model, const set<int>& sources)
			: Geodesic::CXin_Wang(model, sources) {
			for (auto source : sources)
			{
				for (int k = 0; k < model.Neigh(source).size(); ++k)
				{
					m_candidateCoveringFacesForEachAncestor[source].insert(model.Edge(model.Neigh(source)[k].first).indexOfFrontFace);
				}
			}
		}
		void Propagate()
		{
			set<int> tmpDestinations(m_destinations);
			ComputeChildrenOfSource();
			bool fFromQueueOfPseudoSources = UpdateTreeDepthBackWithChoice();

			while (!m_QueueForPseudoSources.empty() || !m_QueueForWindows.empty())
			{
				if (m_QueueForWindows.size() > m_nMaxLenOfWindowQueue)
					m_nMaxLenOfWindowQueue = (int)m_QueueForWindows.size();
				if (m_QueueForPseudoSources.size() > m_nMaxLenOfPseudoSourceQueue)
					m_nMaxLenOfPseudoSourceQueue = m_QueueForPseudoSources.size();
				if (m_QueueForWindows.size() + m_QueueForPseudoSources.size() > m_maxLenOfQueue)
					m_maxLenOfQueue = m_QueueForWindows.size() + m_QueueForPseudoSources.size();
				if (fFromQueueOfPseudoSources) //pseudosource
				{
					int indexOfVert = m_QueueForPseudoSources.top().indexOfVert;
					//cout << "Vert = " <<  indexOfVert << endl;
					//newly added! begin
					//m_destinationVerticesOfAncestors[m_InfoAtVertices[indexOfVert].indexOfAncestor].insert(indexOfVert);
					for (int k = 0; k < model.Neigh(indexOfVert).size(); ++k)
					{
						m_candidateCoveringFacesForEachAncestor[m_InfoAtVertices[indexOfVert].indexOfAncestor].insert(model.Edge(model.Neigh(indexOfVert)[k].first).indexOfFrontFace);
					}
					//newly added! end
					m_QueueForPseudoSources.pop();
					if (m_InfoAtVertices[indexOfVert].disUptodate > m_radius)
						break;
					tmpDestinations.erase(indexOfVert);

					if (!m_destinations.empty() && tmpDestinations.empty())
						return;
					ComputeChildrenOfPseudoSource(indexOfVert);
				}
				else
				{
					QuoteWindow quoteW = m_QueueForWindows.top();
					//newly added! begin
					auto e = quoteW.pWindow->indexOfCurEdge;
					//m_destinationVerticesOfAncestors[quoteW.pWindow->indexOfAncestor].insert(model.Edge(e).indexOfLeftVert);
					//m_destinationVerticesOfAncestors[quoteW.pWindow->indexOfAncestor].insert(model.Edge(e).indexOfRightVert);
					//m_destinationVerticesOfAncestors[quoteW.pWindow->indexOfAncestor].insert(model.Edge(e).indexOfOppositeVert);
					m_candidateCoveringFacesForEachAncestor[quoteW.pWindow->indexOfAncestor].insert(model.Edge(e).indexOfFrontFace);
					//newly added! end
					m_QueueForWindows.pop();
					if (quoteW.disUptodate > m_radius)
						break;
					ComputeChildrenOfWindow(quoteW);
					delete quoteW.pWindow;
				}
				fFromQueueOfPseudoSources = UpdateTreeDepthBackWithChoice();
			}
		}
	};


	//OverPropagation(const Model3D::CRichModel& model, const map<int, double>& sources)
	//	: Geodesic::CXin_Wang(model, sources) {}
	//OverPropagation(const Model3D::CRichModel& model, const set<int>& sources)
	//	: Geodesic::CXin_Wang(model, sources) {}

	vector<vector<tuple<int, double, double, double>>> InferOverPropagatedDistancesForGVD(const Model3D::CRichModel& model, const map<int, double>& sources)
	{
		OverPropagation alg(model, sources);
		alg.Execute();
		vector<vector<tuple<int, double, double, double>>> m_distances(model.GetNumOfFaces());
		for (auto mypair : sources)
		{
			map<int, double> single_sources;
			single_sources[mypair.first] = sources.find(mypair.first)->second;
			set<int> destinations;
			for (auto faceID : alg.m_candidateCoveringFacesForEachAncestor[mypair.first])
			{
				for (int k = 0; k < 3; ++k)
				{
					destinations.insert(model.Face(faceID)[k]);
				}
			}
			Geodesic::CXin_Wang alg_single_source(model, single_sources, destinations);
			alg_single_source.Execute();
			for (auto faceID : alg.m_candidateCoveringFacesForEachAncestor[mypair.first])
			{
				m_distances[faceID].push_back(make_tuple(mypair.first,
					alg_single_source.GetDistanceField()[model.Face(faceID)[0]],
					alg_single_source.GetDistanceField()[model.Face(faceID)[1]],
					alg_single_source.GetDistanceField()[model.Face(faceID)[2]]));
			}
		}
		return m_distances;
	}

	vector<vector<tuple<int, double, double, double>>> InferOverPropagatedDistancesForGVD(const Model3D::CRichModel& model, const set<int>& sources)
	{
		OverPropagation alg(model, sources);
		alg.Execute();
#if 0
		///
		{
			int source = 86;
			vector<Model3D::CBaseModel::CFace> faces;
			for (auto faceID : alg.m_destinationFacesOfAncestors[source])
			{
				faces.push_back(model.Face(faceID));
			}
			Model3D::CBaseModel newModel(model.m_Verts, faces);
			newModel.PreprocessVertsAndFacesIntoBaseModel();
			newModel.SaveObjFile("model1.obj");
		}
		{
			int source = 94;
			vector<Model3D::CBaseModel::CFace> faces;
			for (auto faceID : alg.m_destinationFacesOfAncestors[source])
			{
				faces.push_back(model.Face(faceID));
			}
			Model3D::CBaseModel newModel(model.m_Verts, faces);
			newModel.PreprocessVertsAndFacesIntoBaseModel();
			newModel.SaveObjFile("model2.obj");
		}
		///
#endif
		vector<vector<tuple<int, double, double, double>>> m_distances(model.GetNumOfFaces());
		//int _for_ = 0;
		for (auto source : sources)
		{
			/*cout << ++_for_ << endl;*/
			set<int> single_source;
			single_source.insert(source);
			set<int> destinations;
			for (auto faceID : alg.m_candidateCoveringFacesForEachAncestor[source])
			{
				for (int k = 0; k < 3; ++k)
				{
					destinations.insert(model.Face(faceID)[k]);
				}
			}
			Geodesic::CXin_Wang alg_single_source(model, single_source, destinations);
			alg_single_source.Execute();
			for (auto faceID : alg.m_candidateCoveringFacesForEachAncestor[source])
			{
				m_distances[faceID].push_back(make_tuple(source,
					alg_single_source.GetDistanceField()[model.Face(faceID)[0]],
					alg_single_source.GetDistanceField()[model.Face(faceID)[1]],
					alg_single_source.GetDistanceField()[model.Face(faceID)[2]]));
			}
		}
		return m_distances;
	}

	struct pair_hash
	{
		std::size_t operator () (const std::pair<int, int> const& pair) const
		{
			std::size_t h1 = std::hash<int>()(pair.first);
			std::size_t h2 = std::hash<int>()(pair.second);

			return h1 ^ h2;
		}
		/*template <class T1, class T2>
		std::size_t operator () (std::pair<T1, T2> const& pair) const
		{
			std::size_t h1 = std::hash<T1>()(pair.first);
			std::size_t h2 = std::hash<T2>()(pair.second);

			return h1 ^ h2;
		}*/
	};


	vector<vector<tuple<int, double, double, double>>> InferOverPropagatedDistancesForLRVD(const Model3D::CRichModel& model, const vector<tuple<double, double, double, int>>& sources)
	{

		vector<vector<tuple<int, double, double, double>>> m_distances(model.GetNumOfFaces());


		vector<tuple<int, double, double, double>> bestDistances(model.GetNumOfFaces(), make_tuple(-1, DBL_MAX, DBL_MAX, DBL_MAX));

		auto Average = [](const tuple<int, double, double, double>& mytuple)
		{
			return (get<1>(mytuple) + get<2>(mytuple) + get<3>(mytuple)) / 3.0;
		};

		struct Evt
		{
			int whichSource;
			int toWhichFace;

			double d1, d2, d3;
			const bool operator>(const Evt& other) const
			{
				return d1 + d2 + d3 > other.d1 + other.d2 + other.d3;
			}
			const double GetAverageDistance() const
			{
				return (d1 + d2 + d3) / 3.0;
			}
		};


		//vector<unordered_set<pair<int, int>, pair_hash>> considered(sources.size());;
		vector<tsl::robin_set<int>> considered(sources.size());
		queue<Evt> pending;

		Evt evt;
		for (int i = 0; i < sources.size(); ++i)
		{
			//������ü���Ϳ���
			evt.whichSource = i;
			evt.toWhichFace = get<3>(sources[evt.whichSource]);
			auto sourcePos = Model3D::CPoint3D(get<0>(sources[i]), get<1>(sources[i]), get<2>(sources[i]));

			evt.d1 = (sourcePos - model.Vert(model.Face(evt.toWhichFace)[0])).Len();
			evt.d2 = (sourcePos - model.Vert(model.Face(evt.toWhichFace)[1])).Len();
			evt.d3 = (sourcePos - model.Vert(model.Face(evt.toWhichFace)[2])).Len();


			pending.push(evt);
			//considered.insert(make_pair(evt.whichSource, evt.toWhichFace));
			considered[evt.whichSource].insert(evt.toWhichFace);
		}


		int threeNeighboringFaces[3];
		double t = GetTickCount64();
		Model3D::CPoint3D sourcePos;


		while (!pending.empty())
		{
			Evt& evt = pending.front();

			if (evt.GetAverageDistance() > Average(bestDistances[evt.toWhichFace]))
			{
				if (evt.d1 < get<1>(bestDistances[evt.toWhichFace])
					|| evt.d2 < get<2>(bestDistances[evt.toWhichFace])
					|| evt.d3 < get<3>(bestDistances[evt.toWhichFace]))
				{
					m_distances[evt.toWhichFace].emplace_back(evt.whichSource, evt.d1, evt.d2, evt.d3);
				}
				else {
					pending.pop();
					continue;
				}
			}
			else
			{
				bestDistances[evt.toWhichFace] = make_tuple(evt.whichSource, evt.d1, evt.d2, evt.d3);
				m_distances[evt.toWhichFace].emplace_back(evt.whichSource, evt.d1, evt.d2, evt.d3);
			}

			const int firstEdge = model.GetEdgeIndexFromTwoVertices(model.Face(evt.toWhichFace)[0], model.Face(evt.toWhichFace)[1]);
			threeNeighboringFaces[0] = model.Edge(model.Edge(firstEdge).indexOfReverseEdge).indexOfFrontFace;
			threeNeighboringFaces[1] = model.Edge(model.Edge(firstEdge).indexOfLeftEdge).indexOfFrontFace;
			threeNeighboringFaces[2] = model.Edge(model.Edge(firstEdge).indexOfRightEdge).indexOfFrontFace;

			sourcePos.x = get<0>(sources[evt.whichSource]);
			sourcePos.y = get<1>(sources[evt.whichSource]);
			sourcePos.z = get<2>(sources[evt.whichSource]);


			for (int i = 0; i < 3; ++i)
			{
				if (considered[evt.whichSource].find(threeNeighboringFaces[i]) != considered[evt.whichSource].end() || threeNeighboringFaces[i] == -1)
					continue;
				evt.toWhichFace = threeNeighboringFaces[i];
				evt.d1 = (sourcePos - model.Vert(model.Face(evt.toWhichFace)[0])).Len();
				evt.d2 = (sourcePos - model.Vert(model.Face(evt.toWhichFace)[1])).Len();
				evt.d3 = (sourcePos - model.Vert(model.Face(evt.toWhichFace)[2])).Len();
				pending.push(evt);
				considered[evt.whichSource].insert(evt.toWhichFace);
			}
			pending.pop();
		}

		t = GetTickCount64() - t;
		cerr << "while time: " << t / 1000.0 << "  seconds..." << endl;

		return m_distances;
	}



	//111
	vector<vector<tuple<int, double, double, double>>> InferOverPropagatedDistancesForFastMarching(const Model3D::CRichModel& model, const vector<pair<int, int>>& sources,map<pair<int, int>, double>& mpdis)
	{
		vector<vector<tuple<int, double, double, double>>> m_distances(model.GetNumOfFaces());
		vector<tuple<int, double, double, double>> bestDistances(model.GetNumOfFaces(), make_tuple(-1, DBL_MAX, DBL_MAX, DBL_MAX));
		auto Average = [](tuple<int, double, double, double> mytuple)
		{
			return (get<1>(mytuple) + get<2>(mytuple) + get<3>(mytuple)) / 3.0;
		};
		struct Evt
		{
			int whichSource;
			int toWhichFace;
			double d1, d2, d3;
			bool operator>(const Evt& other) const
			{
				return d1 + d2 + d3 > other.d1 + other.d2 + other.d3;
			}
			double GetAverageDistance() const
			{
				return (d1 + d2 + d3) / 3.0;
			}
		};


		unordered_set<pair<int, int>, pair_hash> considered;
		queue<Evt> pending;
		//initialization


		for (int i = 0; i < sources.size(); ++i)
		{
			Evt evt;
			//������ü���Ϳ���
			evt.whichSource = i;


			evt.toWhichFace = sources[evt.whichSource].second;

			evt.d1 = (model.Vert(sources[evt.whichSource].first) - model.Vert(model.Face(evt.toWhichFace)[0])).Len();
			evt.d2 = (model.Vert(sources[evt.whichSource].first) - model.Vert(model.Face(evt.toWhichFace)[1])).Len();
			evt.d3 = (model.Vert(sources[evt.whichSource].first) - model.Vert(model.Face(evt.toWhichFace)[2])).Len();
			/*evt.d1 = mpdis[make_pair(sources[evt.whichSource].first, model.Face(evt.toWhichFace)[0])];
			evt.d2 = mpdis[make_pair(sources[evt.whichSource].first, model.Face(evt.toWhichFace)[1])];
			evt.d3 = mpdis[make_pair(sources[evt.whichSource].first, model.Face(evt.toWhichFace)[2])];*/

			pending.push(evt);
			considered.insert(make_pair(evt.whichSource, evt.toWhichFace));

		}




		vector<int> threeNeighboringFaces;
		Evt evt;
		threeNeighboringFaces.resize(3);


		while (!pending.empty())
		{
			//cout << "??? " << endl;

			auto top = pending.front();
			pending.pop();
			if (top.GetAverageDistance() > Average(bestDistances[top.toWhichFace]))
			{
				if (top.d1 > get<1>(bestDistances[top.toWhichFace])
					&& top.d2 > get<2>(bestDistances[top.toWhichFace])
					&& top.d3 > get<3>(bestDistances[top.toWhichFace]))
				{
					continue;
				}
				else
				{
					m_distances[top.toWhichFace].push_back(make_tuple(top.whichSource, top.d1, top.d2, top.d3));
				}
			}
			else
			{
				bestDistances[top.toWhichFace] = make_tuple(top.whichSource, top.d1, top.d2, top.d3);
				m_distances[top.toWhichFace].push_back(make_tuple(top.whichSource, top.d1, top.d2, top.d3));
			}

			//cout << "!!!" << endl;
			int firstEdge = model.GetEdgeIndexFromTwoVertices(model.Face(top.toWhichFace)[0], model.Face(top.toWhichFace)[1]);
			threeNeighboringFaces[0] = model.Edge(model.Edge(firstEdge).indexOfReverseEdge).indexOfFrontFace;
			threeNeighboringFaces[1] = model.Edge(model.Edge(firstEdge).indexOfLeftEdge).indexOfFrontFace;
			threeNeighboringFaces[2] = model.Edge(model.Edge(firstEdge).indexOfRightEdge).indexOfFrontFace;


			evt.whichSource = top.whichSource;



			//cout << "..." << endl;
			for (int i = 0; i < 3; ++i)
			{
				if (considered.find(make_pair(top.whichSource, threeNeighboringFaces[i])) != considered.end())
					continue;

				if (threeNeighboringFaces[i] == -1) {
					continue;
				}

				evt.toWhichFace = threeNeighboringFaces[i];
				/*evt.d1 = mpdis[{sources[top.whichSource].first, model.Face(evt.toWhichFace)[0]}];
				evt.d2 = mpdis[{sources[top.whichSource].first, model.Face(evt.toWhichFace)[1]}];
				evt.d3 = mpdis[{sources[top.whichSource].first, model.Face(evt.toWhichFace)[2]}];*/
				evt.d1 = (model.Vert(sources[evt.whichSource].first) - model.Vert(model.Face(evt.toWhichFace)[0])).Len();
				evt.d2 = (model.Vert(sources[evt.whichSource].first) - model.Vert(model.Face(evt.toWhichFace)[1])).Len();
				evt.d3 = (model.Vert(sources[evt.whichSource].first) - model.Vert(model.Face(evt.toWhichFace)[2])).Len();
				pending.push(evt);
				considered.insert(make_pair(evt.whichSource, evt.toWhichFace));
			}
		}
		return m_distances;
	}



	//no use....
	template <class HDS>
	struct Build_prism : public CGAL::Modifier_base<HDS> {
		Point_3 v1, v2, v3;
		Vector_3 normal;
		double height;
		Build_prism() {}
		Build_prism(Point_3 v1, Point_3 v2, Point_3 v3, Vector_3 normal, double h = 10000)
			: v1(v1), v2(v2), v3(v3), normal(normal), height(h) {}
		void operator()(HDS& hds) {
			// Postcondition: hds is a valid polyhedral surface.
			CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
			B.begin_surface(6, 8, 24);
			typedef typename HDS::Vertex   Vertex;
			typedef typename Vertex::Point Point;
			B.add_vertex(v1);
			B.add_vertex(v2);
			B.add_vertex(v3);
			B.add_vertex(v1 + height * normal);
			B.add_vertex(v2 + height * normal);
			B.add_vertex(v3 + height * normal);
			B.begin_facet();
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(1);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(3);
			B.add_vertex_to_facet(4);
			B.add_vertex_to_facet(5);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(1);
			B.add_vertex_to_facet(3);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(3);
			B.add_vertex_to_facet(1);
			B.add_vertex_to_facet(4);
			B.end_facet();
			//
			B.begin_facet();
			B.add_vertex_to_facet(1);
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(4);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(4);
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(5);
			B.end_facet();
			//
			B.begin_facet();
			B.add_vertex_to_facet(2);
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(5);
			B.end_facet();
			B.begin_facet();
			B.add_vertex_to_facet(5);
			B.add_vertex_to_facet(0);
			B.add_vertex_to_facet(3);
			B.end_facet();
			//
			B.end_surface();
		}
	};
	//no use....
	void Test_Build_prism()
	{
		Vector_3 normal(0, 0, 1);
		Point_3 v1(0, 0, 0);
		Point_3 v2(1, 0, 0);
		Point_3 v3(0, 1, 0);
		Polyhedron P;
		Build_prism<HalfedgeDS> prism(v1, v2, v3, normal, 10);
		P.delegate(prism);
		ofstream out("prism.off");
		out << P;
		out.close();
	}

	struct Triangle_Base_Convex_Top
	{
		const Model3D::CRichModel& model;
		double x1, y1;
		double x2, y2;
		double x3, y3;
		int faceID;
		//from 2d to 3d
		Model3D::CPoint3D Get3DPoint(double x, double y) const
		{
			if (faceID == -1)
			{
				return Model3D::CPoint3D(x, y, 0);
			}
			{
				auto coord = model.m_coord_systems_faces[faceID];
				return get<0>(coord) + x * get<1>(coord) + y * get<2>(coord);
			}
		}
		struct MyFacet
		{
			bool isWall;
			int ancestor;
			//double h1, h2, h3;
			int wallID;
			double a, b, c;
			MyFacet()
			{
				isWall = true;
				wallID = -1;
				ancestor = -1;
			}
			MyFacet(const Triangle_Base_Convex_Top* tble, double h1, double h2, double h3, int ancestor)
				: ancestor(ancestor)
			{
				isWall = false;
				Eigen::Matrix3d m;
				m << tble->x1, tble->y1, 1,
					tble->x2, tble->y2, 1,
					tble->x3, tble->y3, 1;
				Eigen::Vector3d right;
				right << h1, h2, h3;
				auto res = m.inverse() * right;
				a = res(0);
				b = res(1);
				c = res(2);
			}
			bool OnUpperSide(double x, double y, double z) const
			{
				return z > a * x + b * y + c;
			}
		};
		struct MyVertex {
			int facet1;
			int facet2;
			int facet3;
			double x, y, h;
			bool IsBottomVertex() const
			{
				return facet3 == -1 || facet1 == -1 || facet2 == -1;
			}
		};

		struct MyEdge
		{
			int vertex1;
			int vertex2;
			int facet1;
			int facet2;
		};

		vector<MyFacet> m_linearFields;
		vector<MyVertex> vertexPool;
		set<int> survivingVertices;
		vector<MyEdge> edgePool;
		set<int> survivingEdges;
		//set<int> survivingFacets;

		struct SubFacet
		{
			int faceID;
			bool isComplete;
			Polygon_2 boundary; //note: (0, 0), (l, 0), (.x, .y)
		};

		Triangle_Base_Convex_Top(double x1, double y1,
			double x2, double y2,
			double x3, double y3) : x1(x1), y1(y1), x2(x2), y2(y2), x3(x3), y3(y3),
			model(Model3D::CRichModel("")), faceID(-1)
		{
			MyFacet facet;
			facet.isWall = true;
			for (int i = 0; i < 3; ++i)
			{
				facet.wallID = i;
				m_linearFields.push_back(facet);
				//survivingFacets.insert(m_linearFields.size() - 1);
			}
			double maxHeight = 100000;//model.m_maxEdgeLength * model.GetNumOfEdges() / 2 + 10000;
			m_linearFields.push_back(MyFacet(this, maxHeight, maxHeight, maxHeight, -1));
			//survivingFacets.insert(m_linearFields.size() - 1);

			MyVertex v;
			v.x = x1;
			v.y = y1;
			v.h = maxHeight;
			v.facet1 = 2;
			v.facet2 = 0;
			v.facet3 = 3;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x2;
			v.y = y2;
			v.h = maxHeight;
			v.facet1 = 0;
			v.facet2 = 1;
			v.facet3 = 3;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x3;
			v.y = y3;
			v.h = maxHeight;
			v.facet1 = 1;
			v.facet2 = 2;
			v.facet3 = 3;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);

			v.x = x1;
			v.y = y1;
			v.h = -maxHeight;
			v.facet1 = 2;
			v.facet2 = 0;
			v.facet3 = -1;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x2;
			v.y = y2;
			v.h = -maxHeight;
			v.facet1 = 0;
			v.facet2 = 1;
			v.facet3 = -1;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x3;
			v.y = y3;
			v.h = -maxHeight;
			v.facet1 = 1;
			v.facet2 = 2;
			v.facet3 = -1;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);

			MyEdge edge;
			edge.vertex1 = 0;
			edge.vertex2 = 1;
			edge.facet1 = 0;
			edge.facet2 = 3;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);
			edge.vertex1 = 1;
			edge.vertex2 = 2;
			edge.facet1 = 1;
			edge.facet2 = 3;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);
			edge.vertex1 = 2;
			edge.vertex2 = 0;
			edge.facet1 = 2;
			edge.facet2 = 3;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);

			edge.vertex1 = 0;
			edge.vertex2 = 3;
			edge.facet1 = 2;
			edge.facet2 = 0;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);


			edge.vertex1 = 1;
			edge.vertex2 = 4;
			edge.facet1 = 0;
			edge.facet2 = 1;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);

			edge.vertex1 = 2;
			edge.vertex2 = 5;
			edge.facet1 = 1;
			edge.facet2 = 2;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);
		}

		Triangle_Base_Convex_Top(const Model3D::CRichModel& model,
			int faceID) : faceID(faceID), model(model)
		{
			x1 = 0; y1 = 0;
			int edgeID = model.GetEdgeIndexFromTwoVertices(model.Face(faceID)[0], model.Face(faceID)[1]);
			x2 = model.Edge(edgeID).length;
			y2 = 0;
			x3 = model.Edge(edgeID).coordOfOppositeVert.x();
			y3 = model.Edge(edgeID).coordOfOppositeVert.y();

			MyFacet facet;
			facet.isWall = true;

			//加入三个wall
			for (int i = 0; i < 3; ++i)
			{
				facet.wallID = i;
				m_linearFields.push_back(facet);
				//survivingFacets.insert(m_linearFields.size() - 1);
			}
			double maxHeight = model.m_maxEdgeLength * model.GetNumOfEdges() / 2 + 10000;

			//给三棱柱封顶，最高的那个平面
			m_linearFields.push_back(MyFacet(this, maxHeight, maxHeight, maxHeight, -1));
			//survivingFacets.insert(m_linearFields.size() - 1);

			MyVertex v;
			v.x = x1;
			v.y = y1;
			v.h = maxHeight;
			v.facet1 = 2;
			v.facet2 = 0;
			v.facet3 = 3;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x2;
			v.y = y2;
			v.h = maxHeight;
			v.facet1 = 0;
			v.facet2 = 1;
			v.facet3 = 3;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x3;
			v.y = y3;
			v.h = maxHeight;
			v.facet1 = 1;
			v.facet2 = 2;
			v.facet3 = 3;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);

			v.x = x1;
			v.y = y1;
			v.h = -maxHeight;
			v.facet1 = 2;
			v.facet2 = 0;
			v.facet3 = -1;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x2;
			v.y = y2;
			v.h = -maxHeight;
			v.facet1 = 0;
			v.facet2 = 1;
			v.facet3 = -1;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);
			v.x = x3;
			v.y = y3;
			v.h = -maxHeight;
			v.facet1 = 1;
			v.facet2 = 2;
			v.facet3 = -1;
			vertexPool.push_back(v);
			survivingVertices.insert(vertexPool.size() - 1);

			MyEdge edge;
			edge.vertex1 = 0;
			edge.vertex2 = 1;
			edge.facet1 = 0;
			edge.facet2 = 3;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);
			edge.vertex1 = 1;
			edge.vertex2 = 2;
			edge.facet1 = 1;
			edge.facet2 = 3;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);
			edge.vertex1 = 2;
			edge.vertex2 = 0;
			edge.facet1 = 2;
			edge.facet2 = 3;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);

			edge.vertex1 = 0;
			edge.vertex2 = 3;
			edge.facet1 = 2;
			edge.facet2 = 0;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);


			edge.vertex1 = 1;
			edge.vertex2 = 4;
			edge.facet1 = 0;
			edge.facet2 = 1;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);

			edge.vertex1 = 2;
			edge.vertex2 = 5;
			edge.facet1 = 1;
			edge.facet2 = 2;
			edgePool.push_back(edge);
			survivingEdges.insert(edgePool.size() - 1);
		}

		void AddFacet(const Triangle_Base_Convex_Top::MyFacet& facet_new)
		{
			set<int> uselessVertexIDs;
			for (auto v : survivingVertices)
			{
				auto v_pos = vertexPool[v];
				if (facet_new.OnUpperSide(v_pos.x, v_pos.y, v_pos.h))
					uselessVertexIDs.insert(v);
			}
			if (uselessVertexIDs.empty())
				return;
			for (auto v : uselessVertexIDs)
				survivingVertices.erase(v);


			//һ�����Կռ�
			m_linearFields.push_back(facet_new);
			int facet_new_id = m_linearFields.size() - 1;
			map<int, set<int>> verticesInFacet; //each site facet contains two new vertices

			set<int> uselessEdges;
			for (auto eID : survivingEdges)
			{
				//��ǰ��������Ҵ�
				bool flag1 = (survivingVertices.find(edgePool[eID].vertex1) != survivingVertices.end());
				bool flag2 = (survivingVertices.find(edgePool[eID].vertex2) != survivingVertices.end());
				if (flag1 && flag2)
					continue;

				//ֻ��һ�����Ҵ�
				if (flag1 && !flag2)
				{

					//һ���µ����������µĶ���
					MyVertex vertex_new;
					vertex_new.facet1 = edgePool[eID].facet1;
					vertex_new.facet2 = edgePool[eID].facet2;
					vertex_new.facet3 = facet_new_id;


					double delta1 = facet_new.a * vertexPool[edgePool[eID].vertex1].x
						+ facet_new.b * vertexPool[edgePool[eID].vertex1].y
						+ facet_new.c - vertexPool[edgePool[eID].vertex1].h;

					double delta2 = facet_new.a * vertexPool[edgePool[eID].vertex2].x
						+ facet_new.b * vertexPool[edgePool[eID].vertex2].y
						+ facet_new.c - vertexPool[edgePool[eID].vertex2].h;



					double lambda = (delta1 - 0) / (delta1 - delta2);
					vertex_new.x = (1 - lambda) * vertexPool[edgePool[eID].vertex1].x
						+ lambda * vertexPool[edgePool[eID].vertex2].x;
					vertex_new.y = (1 - lambda) * vertexPool[edgePool[eID].vertex1].y
						+ lambda * vertexPool[edgePool[eID].vertex2].y;
					vertex_new.h = (1 - lambda) * vertexPool[edgePool[eID].vertex1].h
						+ lambda * vertexPool[edgePool[eID].vertex2].h;
					vertexPool.push_back(vertex_new);


					if (isnan(vertex_new.x)) {
						cout << "AddFacet x nan wrony" << endl;
						exit(-1);
					}
					if (isnan(vertex_new.y)) {
						cout << "AddFacet y nan wrony" << endl;
						exit(-1);
					}
					if (isnan(vertex_new.h)) {
						cout << "AddFacet h nan wrony" << endl;
						exit(-1);
					}


					int vertex_new_id = vertexPool.size() - 1;
					survivingVertices.insert(vertex_new_id);
					edgePool[eID].vertex2 = vertex_new_id;
					verticesInFacet[edgePool[eID].facet1].insert(vertex_new_id);
					verticesInFacet[edgePool[eID].facet2].insert(vertex_new_id);
					//这地方是不是少了一个？
				}
				else if (!flag1 && flag2)
				{
					MyVertex vertex_new;
					vertex_new.facet1 = edgePool[eID].facet1;
					vertex_new.facet2 = edgePool[eID].facet2;
					vertex_new.facet3 = facet_new_id;
					double delta1 = facet_new.a * vertexPool[edgePool[eID].vertex1].x
						+ facet_new.b * vertexPool[edgePool[eID].vertex1].y
						+ facet_new.c - vertexPool[edgePool[eID].vertex1].h;
					double delta2 = facet_new.a * vertexPool[edgePool[eID].vertex2].x
						+ facet_new.b * vertexPool[edgePool[eID].vertex2].y
						+ facet_new.c - vertexPool[edgePool[eID].vertex2].h;

					double lambda = (delta1 - 0) / (delta1 - delta2);
					vertex_new.x = (1 - lambda) * vertexPool[edgePool[eID].vertex1].x
						+ lambda * vertexPool[edgePool[eID].vertex2].x;
					vertex_new.y = (1 - lambda) * vertexPool[edgePool[eID].vertex1].y
						+ lambda * vertexPool[edgePool[eID].vertex2].y;
					vertex_new.h = (1 - lambda) * vertexPool[edgePool[eID].vertex1].h
						+ lambda * vertexPool[edgePool[eID].vertex2].h;


					if (isnan(vertex_new.x)) {
						cout << "AddFacet x nan wrony" << endl;
						exit(-1);
					}
					if (isnan(vertex_new.y)) {
						cout << "AddFacet y nan wrony" << endl;
						exit(-1);
					}
					if (isnan(vertex_new.h)) {
						cout << "AddFacet h nan wrony" << endl;
						exit(-1);
					}


					vertexPool.push_back(vertex_new);
					int vertex_new_id = vertexPool.size() - 1;
					survivingVertices.insert(vertex_new_id);
					edgePool[eID].vertex1 = vertex_new_id;
					verticesInFacet[edgePool[eID].facet1].insert(vertex_new_id);
					verticesInFacet[edgePool[eID].facet2].insert(vertex_new_id);
				}
				else// (!flag1 && !flag2)
				{
					uselessEdges.insert(eID);
				}
			}
			for (auto eID : uselessEdges)
				survivingEdges.erase(eID);


			for (auto mypair : verticesInFacet)
			{
				MyEdge edge_new;
				edge_new.vertex1 = *mypair.second.begin();
				edge_new.vertex2 = *mypair.second.rbegin();
				edge_new.facet2 = facet_new_id;
				if (vertexPool[edge_new.vertex1].facet1 == vertexPool[edge_new.vertex2].facet1)
				{
					edge_new.facet1 = vertexPool[edge_new.vertex1].facet1;
				}
				else if (vertexPool[edge_new.vertex1].facet1 == vertexPool[edge_new.vertex2].facet2)
				{
					edge_new.facet1 = vertexPool[edge_new.vertex1].facet1;
				}
				else if (vertexPool[edge_new.vertex1].facet2 == vertexPool[edge_new.vertex2].facet1)
				{
					edge_new.facet1 = vertexPool[edge_new.vertex1].facet2;
				}
				else //if (vertexPool[edge_new.vertex1].facet2 == vertexPool[edge_new.vertex2].facet2)
				{
					edge_new.facet1 = vertexPool[edge_new.vertex1].facet2;
				}
				edgePool.push_back(edge_new);
				survivingEdges.insert(edgePool.size() - 1);
			}
		}

		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> GetSegments() const
		{
			vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> result;
			for (auto eID : survivingEdges)
			{
				auto e = edgePool[eID];
				if (e.facet1 <= 3 || e.facet2 <= 3) //walls or the initial top face
					continue;
				if (isnan(vertexPool[e.vertex1].x) || isnan(vertexPool[e.vertex1].y) || isnan(vertexPool[e.vertex2].x) || isnan(vertexPool[e.vertex2].y)) {
					cout << "GetSegments function wrong" << endl;
					exit(-1);
				}
				auto p1 = Get3DPoint(vertexPool[e.vertex1].x, vertexPool[e.vertex1].y);
				auto p2 = Get3DPoint(vertexPool[e.vertex2].x, vertexPool[e.vertex2].y);
				result.push_back(make_pair(p1, p2));
			}
			return result;
		}
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> GetSegments(tsl::robin_map <int, int>& mp) const
		{
			vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> result;
			for (auto eID : survivingEdges)
			{
				auto e = edgePool[eID];
				if (e.facet1 <= 3 || e.facet2 <= 3) //walls or the initial top face
					continue;
				if (isnan(vertexPool[e.vertex1].x) || isnan(vertexPool[e.vertex1].y) || isnan(vertexPool[e.vertex2].x) || isnan(vertexPool[e.vertex2].y)) {
					cout << "GetSegments function wrong" << endl;
					exit(-1);
				}
				mp[e.vertex1]++;
				mp[e.vertex2]++;
				auto p1 = Get3DPoint(vertexPool[e.vertex1].x, vertexPool[e.vertex1].y);
				auto p2 = Get3DPoint(vertexPool[e.vertex2].x, vertexPool[e.vertex2].y);
				result.push_back(make_pair(p1, p2));
			}
			return result;
		}

		void RegisterSubFacets(map<int, vector<SubFacet>>& facetsBySites) const
		{
			map<int, vector<Point_2>> verticesBySites;
			for (auto vID : survivingVertices)
			{
				auto vertex = vertexPool[vID];
				if (vertex.IsBottomVertex())
					continue;
				if (!m_linearFields[vertex.facet1].isWall)
				{
					verticesBySites[m_linearFields[vertex.facet1].ancestor].push_back(Point_2(vertex.x, vertex.y));
				}
				if (!m_linearFields[vertex.facet2].isWall)
				{
					verticesBySites[m_linearFields[vertex.facet2].ancestor].push_back(Point_2(vertex.x, vertex.y));
				}
				if (!m_linearFields[vertex.facet3].isWall)
				{
					verticesBySites[m_linearFields[vertex.facet3].ancestor].push_back(Point_2(vertex.x, vertex.y));
				}
			}
			for (auto mypair : verticesBySites)
			{
				int site = mypair.first;
				//optional
				//if (P.area() < 0)
				//	P.reverse_orientation();
				//struct SubFacet
				//{
				//	int faceID;
				//	bool isComplete;
				//	Polygon_2 boundary; //note: (0, 0), (l, 0), (.x, .y)
				//};
				SubFacet facet;
				facet.faceID = faceID;
				if (survivingVertices.size() <= 6)//top: three vertices; bottom: three vertices
				{
					facet.isComplete = true;
				}
				else
				{
					facet.isComplete = false;
				}
				CGAL::convex_hull_2(mypair.second.begin(), mypair.second.end(), back_inserter(facet.boundary));
				facetsBySites[site].push_back(facet);
			}
		}

		Polyhedron ExtractLowerEnvelop() const //only for visualization
		{
			map<Point_3, int> fromCGALVertex2MyID;
			vector<Point_3> vertices;
			for (auto v : survivingVertices)
			{
				Point_3 v_cgal(vertexPool[v].x, vertexPool[v].y, vertexPool[v].h);
				fromCGALVertex2MyID[v_cgal] = v;
				vertices.push_back(v_cgal);
			}
			//vertices.push_back(Point_3(x1, y1, 0));
			//fromCGALVertex2MyID[vertices.back()] = -1;
			//vertices.push_back(Point_3(x2, y2, 0));
			//fromCGALVertex2MyID[vertices.back()] = -1;
			//vertices.push_back(Point_3(x3, y3, 0));
			//fromCGALVertex2MyID[vertices.back()] = -1;
			//convex hull;
			Polyhedron poly;
			CGAL::convex_hull_3(vertices.begin(), vertices.end(), poly);
			return poly;
		}
	};

	void WriteLineObjFile(const vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>>& segs, const char* filename)
	{
		ofstream out(filename);
		int id(0);
		for (int i = 0; i < segs.size(); ++i)
		{
			out << "v " << segs[i].first.x << " " << segs[i].first.y << " " << segs[i].first.z << endl;
			out << "v " << segs[i].second.x << " " << segs[i].second.y << " " << segs[i].second.z << endl;
			out << "l " << id + 1 << " " << id + 2 << endl;
			id += 2;
		}
		out.close();
	}

	map<int, vector<Triangle_Base_Convex_Top::SubFacet>> GetRVD_Regions(const Model3D::CRichModel& model, const vector<tuple<double, double, double, int>>& sources)
	{
		map<int, vector<Triangle_Base_Convex_Top::SubFacet>> regions;
		auto resultingField = InferOverPropagatedDistancesForLRVD(model, sources);
#pragma omp parallel for num_threads(4)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				Triangle_Base_Convex_Top::SubFacet facet;
				facet.faceID = faceID;
				facet.isComplete = true;
#pragma omp critical
				{
					regions[get<0>(resultingField[faceID][0])].push_back(facet);
				}
			}
			else
			{
				Triangle_Base_Convex_Top tble(model, faceID);
				for (int j = 0; j < resultingField[faceID].size(); ++j)
				{
					Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
						get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
						get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
						get<0>(resultingField[faceID][j]));
					tble.AddFacet(facet);
				}
#pragma omp critical
				{
					tble.RegisterSubFacets(regions);
				}
			}
		}
		return regions;
	}



	vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> GetFastMarching_Bisectors(const Model3D::CRichModel& model, const vector<pair<int, int>>& sources, map<pair<int, int>, double>& mpdis)
	{
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> segs;

		double t = GetTickCount64();
		cout << "infer start" << endl;
		cout << "before mpsize " << mpdis.size() << endl;
		auto resultingField = InferOverPropagatedDistancesForFastMarching(model, sources, mpdis);
		cout << "after mpsize " << mpdis.size() << endl;
		t = GetTickCount64() - t;
		cerr << "InferOverPropagatedDistancesForFastMarching time: " << t / 1000.0 << "  seconds..." << endl;

		t = GetTickCount64();
#pragma omp parallel for num_threads(4)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				continue;
			}
			Triangle_Base_Convex_Top tble(model, faceID);
			for (int j = 0; j < resultingField[faceID].size(); ++j)
			{
				Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
					get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
					get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
					get<0>(resultingField[faceID][j]));
				tble.AddFacet(facet);
			}
			auto segs_face = tble.GetSegments();
#pragma omp critical
			{
				copy(segs_face.begin(), segs_face.end(), back_inserter(segs));
			}
		}
		t = GetTickCount64() - t;
		cerr << "Triangle_Base_Convex_Top for all the faces: " << t / 1000.0 << "  seconds..." << endl;
		return segs;
	}

	map<int, vector<Triangle_Base_Convex_Top::SubFacet>> GetGVD_Regions(const Model3D::CRichModel& model, const map<int, double>& sources)
	{
		map<int, vector<Triangle_Base_Convex_Top::SubFacet>> regions;
		auto resultingField = InferOverPropagatedDistancesForGVD(model, sources);
#pragma omp parallel for num_threads(4)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				Triangle_Base_Convex_Top::SubFacet facet;
				facet.faceID = faceID;
				facet.isComplete = true;
#pragma omp critical
				{
					regions[get<0>(resultingField[faceID][0])].push_back(facet);
				}
			}
			else
			{
				Triangle_Base_Convex_Top tble(model, faceID);
				for (int j = 0; j < resultingField[faceID].size(); ++j)
				{
					Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
						get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
						get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
						get<0>(resultingField[faceID][j]));
					tble.AddFacet(facet);
				}
#pragma omp critical
				{
					tble.RegisterSubFacets(regions);
				}
			}
		}
		return regions;
	}




	pair<vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>>,pair<double,double>> GetLRVD_Bisectors(const Model3D::CRichModel& model, const vector<tuple<double, double, double, int>>& sources)
	{
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> segs;

		double t = GetTickCount64();
		cout << "infer start" << endl;
		auto resultingField = InferOverPropagatedDistancesForLRVD(model, sources);
		t = GetTickCount64() - t;
		double infer = t / 1000.0;
		cerr << "InferOverPropagatedDistancesForLRVD time: " << t / 1000.0 << "  seconds..." << endl;

		t = GetTickCount64();


#pragma omp parallel for num_threads(8)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				continue;
			}
			Triangle_Base_Convex_Top tble(model, faceID);
			for (int j = 0; j < resultingField[faceID].size(); ++j)
			{
				Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
					get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
					get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
					get<0>(resultingField[faceID][j]));
				tble.AddFacet(facet);
			}
			auto segs_face = tble.GetSegments();
#pragma omp critical
			{
				copy(segs_face.begin(), segs_face.end(), back_inserter(segs));
			}
		}
		t = GetTickCount64() - t;
		double cutt = t/1000;
		cerr << "Triangle_Base_Convex_Top for all the faces: " << t / 1000.0 << "  seconds..." << endl;
		return { segs, { infer,cutt } };
	}

	vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> GetLRVD_BisectorsWithDui(const Model3D::CRichModel& model, const vector<tuple<double, double, double, int>>& sources)
	{
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> segs;

		double t = GetTickCount64();
		cout << "infer start" << endl;
		auto resultingField = InferOverPropagatedDistancesForLRVD(model, sources);
		t = GetTickCount64() - t;
		cerr << "InferOverPropagatedDistancesForLRVD time: " << t / 1000.0 << "  seconds..." << endl;

		t = GetTickCount64();

		vector<Model3D::CPoint3D> sources3d(sources.size());

		for(int i=0;i<sources.size();++i)
			sources3d[i]= Model3D::CPoint3D(get<0>(sources[i]), get<1>(sources[i]), get<2>(sources[i]));

		auto cross = [](const Model3D::CPoint3D& a, const Model3D::CPoint3D& b) {
			return Model3D::CPoint3D(a.y * b.z - b.y * a.z, -a.x * b.z + a.z * b.x, a.x * b.y - a.y * b.x);
		};

		auto dot = [](const Model3D::CPoint3D& a, const Model3D::CPoint3D& b) {
			return a.x * b.x + a.y * b.y + a.z * b.z;
		};


		ofstream out("other.obj");
		for (int i = 0; i < sources.size(); ++i) {
			out << "v " << get<0>(sources[i]) << " " << get<1>(sources[i]) << " " << get<2>(sources[i]) << endl;
		}

#pragma omp parallel for num_threads(8)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				continue;
			}
			Triangle_Base_Convex_Top tble(model, faceID);
			for (int j = 0; j < resultingField[faceID].size(); ++j)
			{
				Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
					get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
					get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
					get<0>(resultingField[faceID][j]));
				tble.AddFacet(facet);
			}
			tsl::robin_map <int, int>mp;
			auto segs_face = tble.GetSegments(mp);

#pragma omp critical
			{

				for (auto v : mp) {
					if (v.second == 3)
					{
						auto face_id = v.first;
						if (tble.vertexPool[face_id].facet3 != -1 && tble.vertexPool[face_id].facet1 != -1 && tble.vertexPool[face_id].facet2 != -1) {
							int a = tble.m_linearFields[tble.vertexPool[face_id].facet1].ancestor;
							int b = tble.m_linearFields[tble.vertexPool[face_id].facet2].ancestor;
							int c = tble.m_linearFields[tble.vertexPool[face_id].facet3].ancestor;

							//��ʼ������face���ж���
							auto n1 = cross(model.Vert(model.Face(faceID)[0]) - model.Vert(model.Face(faceID)[1]), model.Vert(model.Face(faceID)[0]) - model.Vert(model.Face(faceID)[2]));
							auto n2 = cross(sources3d[a] - sources3d[b], sources3d[a] - sources3d[c]);
							if (dot(n2, n1) > 0) {
								out << "f " << a + 1 << " " << b + 1 << " " << c + 1 << endl;
							}
							else out << "f " << c + 1 << " " << b + 1 << " " << a + 1 << endl;
						}
					}

				}
				copy(segs_face.begin(), segs_face.end(), back_inserter(segs));
			}
		}
		t = GetTickCount64() - t;
		cerr << "Triangle_Base_Convex_Top for all the faces: " << t / 1000.0 << "  seconds..." << endl;
		return segs;
	}


	vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> GetGVD_Bisectors(const Model3D::CRichModel& model, const map<int, double>& sources)
	{
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> segs;
		auto resultingField = InferOverPropagatedDistancesForGVD(model, sources);
#pragma omp parallel for num_threads(4)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				continue;
			}
			Triangle_Base_Convex_Top tble(model, faceID);
			for (int j = 0; j < resultingField[faceID].size(); ++j)
			{
				Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
					get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
					get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
					get<0>(resultingField[faceID][j]));
				tble.AddFacet(facet);
			}
			auto segs_face = tble.GetSegments();
#pragma omp critical
			{
				copy(segs_face.begin(), segs_face.end(), back_inserter(segs));
			}
		}
		return segs;
	}

	map<int, vector<Triangle_Base_Convex_Top::SubFacet>> GetGVD_Regions(const Model3D::CRichModel& model, const set<int>& sources)
	{
		map<int, vector<Triangle_Base_Convex_Top::SubFacet>> regions;
		auto resultingField = InferOverPropagatedDistancesForGVD(model, sources);
#pragma omp parallel for num_threads(4)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				Triangle_Base_Convex_Top::SubFacet facet;
				facet.faceID = faceID;
				facet.isComplete = true;
#pragma omp critical
				{
					regions[get<0>(resultingField[faceID][0])].push_back(facet);
				}
			}
			else
			{
				Triangle_Base_Convex_Top tble(model, faceID);
				for (int j = 0; j < resultingField[faceID].size(); ++j)
				{
					Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
						get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
						get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
						get<0>(resultingField[faceID][j]));
					tble.AddFacet(facet);
				}
#pragma omp critical
				{
					tble.RegisterSubFacets(regions);
				}
			}
		}
		return regions;
	}

	vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> GetGVD_Bisectors(const Model3D::CRichModel& model, const set<int>& sources)
	{
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> segs;
		double t = GetTickCount64();
		auto resultingField = InferOverPropagatedDistancesForGVD(model, sources);
		t = GetTickCount64() - t;
		cerr << "InferOverPropagatedDistancesForGVD time: " << t / 1000.0 << "  seconds..." << endl;

		t = GetTickCount64();
		auto cmp = [](const tuple<int, double, double, double>& a, const tuple<int, double, double, double>& b) {
			return max(get<1>(a), max(get<2>(a), get<3>(a))) < max(get<1>(b), max(get<2>(b), get<3>(b)));
		};

#pragma omp parallel for num_threads(4)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				continue;
			}
			//sort(resultingField[faceID].begin(), resultingField[faceID].end(), cmp);
			Triangle_Base_Convex_Top tble(model, faceID);
			for (int j = 0; j < resultingField[faceID].size(); ++j)
			{
				Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
					get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
					get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
					get<0>(resultingField[faceID][j]));
				tble.AddFacet(facet);
			}
			auto segs_face = tble.GetSegments();
#pragma omp critical
			{
				copy(segs_face.begin(), segs_face.end(), back_inserter(segs));
			}
		}
		t = GetTickCount64() - t;
		cerr << "Triangle_Base_Convex_Top for all faces: " << t / 1000.0 << "  seconds..." << endl;
		return segs;
	}



	//
	//
	//sources get<0>(sources[i]): i�����ڵ����ߵ��±ꣻget<1>(sources[i]): i�����ڵ�face�±�
	vector<vector<tuple<int, double, double, double>>> InferOverPropagatedDistancesForLRVDSeveralPoints(const Model3D::CRichModel& model, const vector<tuple<int, int>>& sources, vector<vector<tuple<double, double, double>>> SSources)
	{

		vector<vector<tuple<int, double, double, double>>> m_distances(model.GetNumOfFaces());


		vector<tuple<int, double, double, double>> bestDistances(model.GetNumOfFaces(), make_tuple(-1, DBL_MAX, DBL_MAX, DBL_MAX));

		auto Average = [](const tuple<int, double, double, double>& mytuple)
		{
			return (get<1>(mytuple) + get<2>(mytuple) + get<3>(mytuple)) / 3.0;
		};

		struct Evt
		{
			int whichSource;
			int toWhichFace;

			double d1, d2, d3;
			const bool operator>(const Evt& other) const
			{
				return d1 + d2 + d3 > other.d1 + other.d2 + other.d3;
			}
			const double GetAverageDistance() const
			{
				return (d1 + d2 + d3) / 3.0;
			}

		};


		auto GetSourcesToPointDis = [&](int sourcesInx, const Model3D::CPoint3D& Point) {
			double dis = 1e18;
			for (auto p : SSources[sourcesInx]) {
				double x = get<0>(p), y = get<1>(p), z = get<2>(p);
				dis = min(dis, sqrt((x - Point.x) * (x - Point.x) + (z - Point.z) * (z - Point.z) + (y - Point.y) * (y - Point.y)));
			}
			return dis;
		};


		vector<tsl::robin_set<int>> considered(sources.size());
		queue<Evt> pending;

		Evt evt;
		//for (int i = 0; i < sources.size(); ++i)
		//{
		//	//������ü���Ϳ���
		//	evt.whichSource = get<0>(sources[i]);
		//	evt.toWhichFace = get<1>(sources[evt.whichSource]);
		//	evt.d1 = GetSourcesToPointDis(evt.whichSource, model.Vert(model.Face(evt.toWhichFace)[0]));
		//	evt.d2 = GetSourcesToPointDis(evt.whichSource, model.Vert(model.Face(evt.toWhichFace)[1]));
		//	evt.d3 = GetSourcesToPointDis(evt.whichSource, model.Vert(model.Face(evt.toWhichFace)[2]));
		//	cout <<"dis: "<< evt.d1 << " " << evt.d2 << " " << evt.d3 << endl;
		//	pending.push(evt);
		//	considered[evt.whichSource].insert(evt.toWhichFace);
		//}


		int threeNeighboringFaces[3];
		double t = GetTickCount64();
		for (int f = 0; f < model.GetNumOfFaces(); ++f)
			for (int i = 0; i < sources.size(); ++i) {
				evt.d1 = GetSourcesToPointDis(i, model.Vert(model.Face(f)[0]));
				evt.d2 = GetSourcesToPointDis(i, model.Vert(model.Face(f)[1]));
				evt.d3 = GetSourcesToPointDis(i, model.Vert(model.Face(f)[2]));
				m_distances[f].emplace_back(i, evt.d1, evt.d2, evt.d3);
			}


		//while (!pending.empty())
		//{
		//	Evt& evt = pending.front();

		//	/*if (evt.whichSource == -1) {
		//		cout << "??? " << endl;
		//	}*/

		//	if (evt.GetAverageDistance() > Average(bestDistances[evt.toWhichFace]))
		//	{
		//		if (evt.d1 < get<1>(bestDistances[evt.toWhichFace])
		//			|| evt.d2 < get<2>(bestDistances[evt.toWhichFace])
		//			|| evt.d3 < get<3>(bestDistances[evt.toWhichFace]))
		//		{
		//			m_distances[evt.toWhichFace].emplace_back(evt.whichSource, evt.d1, evt.d2, evt.d3);
		//		}
		//		else {
		//			pending.pop();
		//			continue;
		//		}
		//	}
		//	else
		//	{
		//		bestDistances[evt.toWhichFace] = make_tuple(evt.whichSource, evt.d1, evt.d2, evt.d3);
		//		m_distances[evt.toWhichFace].emplace_back(evt.whichSource, evt.d1, evt.d2, evt.d3);
		//	}

		//	const int firstEdge = model.GetEdgeIndexFromTwoVertices(model.Face(evt.toWhichFace)[0], model.Face(evt.toWhichFace)[1]);
		//	threeNeighboringFaces[0] = model.Edge(model.Edge(firstEdge).indexOfReverseEdge).indexOfFrontFace;
		//	threeNeighboringFaces[1] = model.Edge(model.Edge(firstEdge).indexOfLeftEdge).indexOfFrontFace;
		//	threeNeighboringFaces[2] = model.Edge(model.Edge(firstEdge).indexOfRightEdge).indexOfFrontFace;


		//	for (int i = 0; i < 3; ++i)
		//	{
		//		if (considered[evt.whichSource].find(threeNeighboringFaces[i]) != considered[evt.whichSource].end() || threeNeighboringFaces[i] == -1)
		//			continue;

		//		evt.toWhichFace = threeNeighboringFaces[i];

		//		evt.d1 = GetSourcesToPointDis(evt.whichSource, model.Vert(model.Face(evt.toWhichFace)[0]));
		//		evt.d2 = GetSourcesToPointDis(evt.whichSource, model.Vert(model.Face(evt.toWhichFace)[1]));
		//		evt.d3 = GetSourcesToPointDis(evt.whichSource, model.Vert(model.Face(evt.toWhichFace)[2]));

		//		pending.push(evt);
		//		considered[evt.whichSource].insert(evt.toWhichFace);
		//	}
		//	pending.pop();
		//}

		t = GetTickCount64() - t;
		cerr << "while time: " << t / 1000.0 << "  seconds..." << endl;

		cout << m_distances.size() << endl;
		/*	for (auto& a : m_distances) {
				cout << a.size() << endl;
			}*/
		return m_distances;
	}

	//
	pair<vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>>,vector<double>> GetMultiSites_Bisectors(const Model3D::CRichModel& model, const vector<tuple<int, int>>& sources, const vector<vector<tuple<double, double, double>>>& SSources)
	{
		vector<pair<Model3D::CPoint3D, Model3D::CPoint3D>> segs;

		double t = GetTickCount64();
		cout << "infer start" << endl;


		auto resultingField = InferOverPropagatedDistancesForLRVDSeveralPoints(model, sources, SSources);
		vector<double> dis(model.GetNumOfVerts(), 1e9);
		for (int fi = 0; fi < model.GetNumOfFaces(); ++fi) {
			//vector<tuple<int, double, double, double>> bestDistances(model.GetNumOfFaces(), make_tuple(-1, DBL_MAX, DBL_MAX, DBL_MAX));
			for (auto d : resultingField[fi])
			{
				//vector<vector<tuple<int, double, double, double>>>
				dis[model.Face(fi)[0]] = min(dis[model.Face(fi)[0]], get<1>(d));
				dis[model.Face(fi)[1]] = min(dis[model.Face(fi)[1]], get<2>(d));
				dis[model.Face(fi)[2]] = min(dis[model.Face(fi)[2]], get<3>(d));
			}	
		}


		t = GetTickCount64() - t;
		cerr << "InferOverPropagatedDistancesForLRVD time: " << t / 1000.0 << "  seconds..." << endl;

		t = GetTickCount64();

#pragma omp parallel for num_threads(8)
		for (int faceID = 0; faceID < model.GetNumOfFaces(); ++faceID)
		{
			if (resultingField[faceID].size() <= 1)
			{
				continue;
			}
			Triangle_Base_Convex_Top tble(model, faceID);
			for (int j = 0; j < resultingField[faceID].size(); ++j)
			{
				Triangle_Base_Convex_Top::MyFacet facet(&tble, get<1>(resultingField[faceID][j]) * get<1>(resultingField[faceID][j]),
					get<2>(resultingField[faceID][j]) * get<2>(resultingField[faceID][j]),
					get<3>(resultingField[faceID][j]) * get<3>(resultingField[faceID][j]),
					get<0>(resultingField[faceID][j]));
				tble.AddFacet(facet);
			}
			auto segs_face = tble.GetSegments();
#pragma omp critical
			{
				copy(segs_face.begin(), segs_face.end(), back_inserter(segs));
			}
		}
		t = GetTickCount64() - t;
		cerr << "Triangle_Base_Convex_Top for all the faces: " << t / 1000.0 << "  seconds..." << endl;
		return { segs,dis };
	}
}

