#pragma once
#include "DistanceApproach.h"
using namespace std;
namespace Geodesic
{
	using namespace Model3D;
	class CDijkstra : public CDistanceApproach
	{
	protected:
		struct InfoAtVertex
		{
			int parent;
			int ancestor;
			int levelOnSequenceTree;
			double disUptodate;
			InfoAtVertex()
			{
				levelOnSequenceTree = -1;
				parent = -1;
				ancestor = -1;
				disUptodate = DBL_MAX;
			}
			InfoAtVertex(double disUptodate, int parent, int ancestor, int level) : disUptodate(disUptodate),
				parent(parent), ancestor(ancestor), levelOnSequenceTree(level)
			{
			}
		};

	protected:
		virtual void Initialize();
		virtual void Dispose();
		virtual void Propagate();
		virtual void CollectExperimentalResults();
	public:
		vector<InfoAtVertex> m_finalResults;
		CDijkstra(const CRichModel& model, int source);
		CDijkstra(const CRichModel& model, const set<int>& sources);
		CDijkstra(const CRichModel& model, const set<int>& sources, const set<int>& destinations);
		CDijkstra(const CRichModel& model, const map<int, double>& sources);
		CDijkstra(const CRichModel& model, const map<int, double>& sources, const set<int>& destinations);
		virtual vector<EdgePoint> BacktraceShortestPath(int end) const;
		virtual int GetAncestor(int vIndex) const;
	};

}