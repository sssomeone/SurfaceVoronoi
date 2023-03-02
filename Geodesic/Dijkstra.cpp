#include "Dijkstra.h"
#include <tuple>
#include <queue>
#include <set>
#include <map>
#include <functional>
using namespace std;
namespace Geodesic
{
	using namespace Model3D;
	CDijkstra::CDijkstra(const CRichModel& model, int source) : CDistanceApproach(model, source)
	{
		m_nameOfAlgorithm = "Dijkstra";
	}

	CDijkstra::CDijkstra(const CRichModel& model, const set<int>& sources) : CDistanceApproach(model, sources)
	{
		m_nameOfAlgorithm = "Dijkstra";
	}

	CDijkstra::CDijkstra(const CRichModel& model, const set<int>& sources, const set<int>& destinations) : CDistanceApproach(model, sources, destinations)
	{
		m_nameOfAlgorithm = "Dijkstra";
	}

	CDijkstra::CDijkstra(const CRichModel& model, const map<int, double>& sources) : CDistanceApproach(model, sources)
	{
		m_nameOfAlgorithm = "Dijkstra";
	}

	CDijkstra::CDijkstra(const CRichModel& model, const map<int, double>& sources, const set<int>& destinations) : CDistanceApproach(model, sources, destinations)
	{
		m_nameOfAlgorithm = "Dijkstra";
	}

	void CDijkstra::Initialize()
	{
		CDistanceApproach::Initialize();
		m_finalResults.resize(model.GetNumOfVerts());
	}

	void CDijkstra::Dispose()
	{
		//nothing...
	}

	void CDijkstra::CollectExperimentalResults()
	{
		m_memory = ((double)model.GetNumOfVerts() * sizeof InfoAtVertex
			+ (double)m_maxLenOfQueue * sizeof(InfoAtVertex*)) / 1024 / 1024;
		for (int i = 0; i < m_finalResults.size(); ++i)
		{
			m_scalarField[i] = m_finalResults[i].disUptodate;
		}
		CDistanceApproach::CollectExperimentalResults();
	}

	void CDijkstra::Propagate()
	{
		set<int> destinations(m_destinations);
		struct Event : std::tuple<double, int, int, int, int>
		{
			Event(double dis, int self, int parent, int root, int level) : std::tuple<double, int, int, int, int>(dis, self, parent, root, level)
			{
			}
			bool operator <(const Event& right) const
			{
				return get<0>(*this) < get<0>(right);
			}
			bool operator >(const Event& right) const
			{
				return get<0>(*this) > get<0>(right);
			}
			double getDistance() const { return get<0>(*this); }
			int getMyID() const { return get<1>(*this); }
			int getDirectParent() const { return get<2>(*this); }
			int getAncestor() const { return get<3>(*this); }
			int getLevel() const { return get<4>(*this); }
		};
		priority_queue<Event, vector<Event>, greater<Event> > evtQue;

		for (map<int, double>::const_iterator it = m_sources.begin();
			it != m_sources.end(); ++it)
		{
			Event evt(it->second, it->first, -1, it->first, 0);
			evtQue.push(evt);
		}
		set<int> fixedSet;

		while (!evtQue.empty())
		{
			if (evtQue.size() > m_maxLenOfQueue)
				m_maxLenOfQueue = evtQue.size();
			Event topEvt = evtQue.top();
			evtQue.pop();
			if (topEvt.getLevel() > m_depthOfResultingTree)
				m_depthOfResultingTree = topEvt.getLevel();
			if (fixedSet.find(topEvt.getMyID()) != fixedSet.end())
				continue;
			fixedSet.insert(topEvt.getMyID());
			m_finalResults[topEvt.getMyID()] = InfoAtVertex(topEvt.getDistance(), topEvt.getDirectParent(), topEvt.getAncestor(), topEvt.getLevel());

			destinations.erase(topEvt.getMyID());
			if (destinations.empty() && !m_destinations.empty())
				break;

			for (int i = 0; i < model.Neigh(topEvt.getMyID()).size(); ++i)
			{
				int v = model.Edge(model.Neigh(topEvt.getMyID())[i].first).indexOfRightVert;
				if (fixedSet.find(v) != fixedSet.end())
					continue;
				double len = model.Edge(model.Neigh(topEvt.getMyID())[i].first).length;
				Event evt(topEvt.getDistance() + len, v, topEvt.getMyID(), topEvt.getAncestor(), topEvt.getLevel() + 1);
				evtQue.push(evt);
			}
		}
	}

	vector<EdgePoint> CDijkstra::BacktraceShortestPath(int end) const
	{
		vector<EdgePoint> path;
		int vIndex = end;
		path.push_back(EdgePoint(end));
		int previousVertex;
		while ((previousVertex = m_finalResults[vIndex].parent) != -1)
		{
			path.push_back(EdgePoint(previousVertex));
			vIndex = previousVertex;
		}
		reverse(path.begin(), path.end());
		return path;
	}

	int CDijkstra::GetAncestor(int vIndex) const
	{
		return m_finalResults[vIndex].ancestor;
	}
}