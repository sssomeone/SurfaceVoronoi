#pragma once
#include "DistanceApproach.h"
namespace Geodesic
{
	using namespace Model3D;
	class CFMM : public CDistanceApproach
	{
	protected:
		struct InfoAtVertex
		{
			bool parentIsVertex;
			int parentIndex;
			int ancestor;
			int levelOnSequenceTree;
			double prop;
			double disUptodate;
			InfoAtVertex()
			{
				parentIsVertex = true;
				parentIndex = -1;
				ancestor = -1;
				levelOnSequenceTree = -1;
				disUptodate = DBL_MAX;
			}
			InfoAtVertex(double disUptodate, bool parentIsVertex, int parentIndex, double prop, int ancestor, int level) : disUptodate(disUptodate),
				parentIndex(parentIndex), ancestor(ancestor), prop(prop), parentIsVertex(parentIsVertex), levelOnSequenceTree(level)
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
		CFMM(const CRichModel& model, int source);
		CFMM(const CRichModel& model, const set<int>& sources);
		CFMM(const CRichModel& model, const set<int>& sources, const set<int>& destinations);
		CFMM(const CRichModel& model, const map<int, double>& sources);
		CFMM(const CRichModel& model, const map<int, double>& sources, const set<int>& destinations);
		virtual vector<EdgePoint> BacktraceShortestPath(int end) const;
		virtual int GetAncestor(int vIndex) const;
	};
}