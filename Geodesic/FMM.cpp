#include "FMM.h"
#include <queue>
#include <tuple>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <functional>
using namespace std;
namespace Geodesic
{
	using namespace Model3D;
	CFMM::CFMM(const CRichModel& model, int source) : CDistanceApproach(model, source)
	{
		m_nameOfAlgorithm = "FMM";
	}

	CFMM::CFMM(const CRichModel& model, const set<int>& sources) : CDistanceApproach(model, sources)
	{
		m_nameOfAlgorithm = "FMM";
	}

	CFMM::CFMM(const CRichModel& model, const set<int>& sources, const set<int>& destinations) : CDistanceApproach(model, sources, destinations)
	{
		m_nameOfAlgorithm = "FMM";
	}

	CFMM::CFMM(const CRichModel& model, const map<int, double>& sources) : CDistanceApproach(model, sources)
	{
		m_nameOfAlgorithm = "FMM";
	}

	CFMM::CFMM(const CRichModel& model, const map<int, double>& sources, const set<int>& destinations) : CDistanceApproach(model, sources, destinations)
	{
		m_nameOfAlgorithm = "FMM";
	}

	void CFMM::Initialize()
	{
		CDistanceApproach::Initialize();
		m_finalResults.resize(model.GetNumOfVerts());
	}

	void CFMM::Dispose()
	{
		//nothing...
	}

	void CFMM::CollectExperimentalResults()
	{
		m_memory = ((double)model.GetNumOfVerts() * sizeof InfoAtVertex
			+ (double)m_maxLenOfQueue * sizeof(InfoAtVertex*)) / 1024 / 1024;
		for (int i = 0; i < m_finalResults.size(); ++i)
		{
			m_scalarField[i] = m_finalResults[i].disUptodate;
		}
		CDistanceApproach::CollectExperimentalResults();
	}

	void CFMM::Propagate()
	{
		set<int> destinations(m_destinations);

		struct Event : std::tuple<double, int, bool, int, double, int, int>
		{
			Event(double dis, int self, bool isVertex, int parent, double prop, int ancestor, int level) : std::tuple<double, int, bool, int, double, int, int>(dis, self, isVertex, parent, prop, ancestor, level)
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
			bool getIsVertex() const { return get<2>(*this); }
			int getParentIndex() const { return get<3>(*this); }
			double getProportion() const { return get<4>(*this); }
			int getAncestor() const { return get<5>(*this); }
			int getLevel() const { return get<6>(*this); }
		};
		priority_queue<Event, vector<Event>, greater<Event> > evtQue;

		for (map<int, double>::const_iterator it = m_sources.begin(); it != m_sources.end(); ++it)
		{
			Event evt(it->second, it->first, true, -1, 0, it->first, 0);
			evtQue.push(evt);
		}

		set<int> fixedVertices;
		while (!evtQue.empty())
		{
			if (evtQue.size() > m_maxLenOfQueue)
				m_maxLenOfQueue = evtQue.size();
			Event topEvt = evtQue.top();
			evtQue.pop();
			if (topEvt.getLevel() > m_depthOfResultingTree)
				m_depthOfResultingTree = topEvt.getLevel();
			int topID = topEvt.getMyID();
			if (fixedVertices.find(topID) != fixedVertices.end())
				continue;
			m_scalarField[topID] = topEvt.getDistance();
			fixedVertices.insert(topID);
			m_finalResults[topID] = InfoAtVertex(topEvt.getDistance(),
				topEvt.getIsVertex(), topEvt.getParentIndex(), topEvt.getProportion(), topEvt.getAncestor(), topEvt.getLevel());
			for (int i = 0; i < (int)model.Neigh(topID).size(); ++i)
			{
				int neigh = model.Edge(model.Neigh(topID)[i].first).indexOfRightVert;
				if (fixedVertices.find(neigh) != fixedVertices.end())
					continue;
				double dis = topEvt.getDistance() + model.Edge(model.Neigh(topID)[i].first).length;
				if (dis < m_scalarField[neigh])
				{
					Event evt(dis, neigh, true, topID, 0, topEvt.getAncestor(), topEvt.getLevel() + 1);
					evtQue.push(evt);
				}
			}
			for (int i = 0; i < (int)model.Neigh(topID).size(); ++i)
			{
				int neigh = model.Edge(model.Neigh(topID)[i].first).indexOfRightVert;
				if (fixedVertices.find(neigh) == fixedVertices.end())
					continue;
				int edge1 = model.Neigh(topID)[i].first;

				if (model.Edge(edge1).indexOfFrontFace != -1)
				{
					int oppoV = model.Edge(edge1).indexOfOppositeVert;
					double angle1 = asin((m_scalarField[neigh] - m_scalarField[topID]) / model.Edge(edge1).length);
					double angle2 = model.Neigh(topID)[i].second;

					double shadow = model.Edge(model.Edge(edge1).indexOfLeftEdge).length * cos(angle1 + angle2);
					if (shadow < 0)
					{
						double dis = m_scalarField[topID] + model.Edge(model.Edge(edge1).indexOfLeftEdge).length;
						if (dis < m_scalarField[oppoV])
						{
							Event evt(dis, oppoV, true, topID, 0, topEvt.getAncestor(), topEvt.getLevel() + 1);
							evtQue.push(evt);
						}
					}
					else if (shadow > model.Edge(edge1).length * cos(angle1))
					{
						double dis = m_scalarField[neigh] + model.Edge(model.Edge(edge1).indexOfRightEdge).length;
						if (dis < m_scalarField[oppoV])
						{
							Event evt(dis, oppoV, true, neigh, 0, topEvt.getAncestor(), topEvt.getLevel() + 1);
							evtQue.push(evt);
						}
					}
					else
					{
						double prop = shadow / (model.Edge(edge1).length * cos(angle1));
						double dis = m_scalarField[topID] + model.Edge(model.Edge(edge1).indexOfLeftEdge).length * sin(angle1 + angle2);
						if (dis < max(m_scalarField[topID], m_scalarField[neigh]))
						{
							dis = max(m_scalarField[topID], m_scalarField[neigh]);
						}
						if (dis < m_scalarField[oppoV])
						{
							Event evt(dis, oppoV, false, edge1, prop, topEvt.getAncestor(), topEvt.getLevel() + 1);
							evtQue.push(evt);
						}
					}
				}
				int edge2 = model.Edge(edge1).indexOfReverseEdge;
				if (model.Edge(edge2).indexOfFrontFace != -1)
				{
					int oppoV = model.Edge(edge2).indexOfOppositeVert;
					double angle2 = model.Neigh(neigh)[model.GetSubindexToVert(neigh, topID)].second;
					double angle1 = asin((m_scalarField[topID] - m_scalarField[neigh]) / model.Edge(edge2).length);

					double shadow = model.Edge(model.Edge(edge2).indexOfLeftEdge).length * cos(angle1 + angle2);
					if (shadow < 0)
					{
						double dis = m_scalarField[neigh] + model.Edge(model.Edge(edge2).indexOfLeftEdge).length;
						if (dis < m_scalarField[oppoV])
						{
							Event evt(dis, oppoV, true, neigh, 0, topEvt.getAncestor(), topEvt.getLevel() + 1);
							evtQue.push(evt);
						}
					}
					else if (shadow > model.Edge(edge2).length * cos(angle1))
					{
						double dis = m_scalarField[topID] + model.Edge(model.Edge(edge2).indexOfRightEdge).length;
						if (dis < m_scalarField[oppoV])
						{
							Event evt(dis, oppoV, true, topID, 0, topEvt.getAncestor(), topEvt.getLevel() + 1);
							evtQue.push(evt);
						}
					}
					else
					{
						double prop = shadow / (model.Edge(edge2).length * cos(angle1));
						double dis = m_scalarField[neigh] + model.Edge(model.Edge(edge2).indexOfLeftEdge).length * sin(angle1 + angle2);
						if (dis < max(m_scalarField[topID], m_scalarField[neigh]))
						{
							dis = max(m_scalarField[topID], m_scalarField[neigh]);
						}
						if (dis < m_scalarField[oppoV])
						{
							Event evt(dis, oppoV, false, edge2, prop, topEvt.getAncestor(), topEvt.getLevel() + 1);
							evtQue.push(evt);
						}
					}
				}
			}
		}
	}

	vector<EdgePoint> CFMM::BacktraceShortestPath(int end) const
	{
		vector<EdgePoint> path;
		EdgePoint ep(end);

		do
		{
			path.push_back(ep);
			EdgePoint parent;
			if (ep.isVertex)
			{
				parent.isVertex = m_finalResults[ep.index].parentIsVertex;
				parent.index = m_finalResults[ep.index].parentIndex;
				parent.proportion = m_finalResults[ep.index].prop;
			}
			else
			{
				int reverseEdgeIndex = model.Edge(ep.index).indexOfReverseEdge;
				assert(reverseEdgeIndex != -1);
				int leftV = model.Edge(reverseEdgeIndex).indexOfLeftVert;
				int rightV = model.Edge(reverseEdgeIndex).indexOfRightVert;
				int oppoV = model.Edge(reverseEdgeIndex).indexOfOppositeVert;
				if (m_finalResults[leftV].parentIsVertex
					&& m_finalResults[rightV].parentIsVertex
					&& m_finalResults[leftV].parentIndex == m_finalResults[rightV].parentIndex)
				{
					assert(m_finalResults[leftV].parentIndex == oppoV);
					parent.isVertex = true;
					parent.index = oppoV;
				}
				else
				{
					if (m_finalResults[leftV].parentIsVertex
						&& m_finalResults[leftV].parentIndex == rightV)
					{
						parent.isVertex = true;
						parent.index = rightV;
					}
					else if (m_finalResults[rightV].parentIsVertex
						&& m_finalResults[rightV].parentIndex == leftV)
					{
						parent.isVertex = true;
						parent.index = leftV;
					}
					else
					{
						//d = ax + by + c;
						//so we are computing a, b, c
						double c = m_finalResults[leftV].disUptodate;
						double a = (m_finalResults[rightV].disUptodate - c) / model.Edge(reverseEdgeIndex).length;
						double b = (m_finalResults[oppoV].disUptodate - c - a * model.Edge(reverseEdgeIndex).coordOfOppositeVert(0)) / model.Edge(reverseEdgeIndex).coordOfOppositeVert(1);
						//gradient: (a, b);
						//starting:
						double scale = model.Edge(reverseEdgeIndex).length
							+ model.Edge(model.Edge(reverseEdgeIndex).indexOfLeftEdge).length;
						scale *= 3;
						Eigen::Vector2d start((1 - ep.proportion) * model.Edge(reverseEdgeIndex).length, 0);
						Eigen::Vector2d end(start(0) - a * scale, start(1) - b * scale);
						double prop = model.ProportionOnLeftEdgeByImage(reverseEdgeIndex,
							end, 1 - ep.proportion);
						if (prop > 0 && prop < 1 - 1e-3)
						{
							parent.isVertex = false;
							parent.index = model.Edge(model.Edge(reverseEdgeIndex).indexOfLeftEdge).indexOfReverseEdge;
							parent.proportion = 1 - prop;
						}
						else if (prop >= 1 - 1e-3 && prop < 1 + 1e-3)
						{
							parent.isVertex = true;
							parent.index = oppoV;
						}
						else
						{
							parent.isVertex = false;
							parent.index = model.Edge(model.Edge(reverseEdgeIndex).indexOfRightEdge).indexOfReverseEdge;
							parent.proportion = 1 - model.ProportionOnRightEdgeByImage(reverseEdgeIndex,
								end, 1 - ep.proportion);
						}
					}
				}
			}
			ep = parent;
		} while (!(ep.isVertex && ep.index == -1));
		reverse(path.begin(), path.end());
		return path;
	}

	int CFMM::GetAncestor(int vIndex) const
	{
		return m_finalResults[vIndex].ancestor;
	}
}