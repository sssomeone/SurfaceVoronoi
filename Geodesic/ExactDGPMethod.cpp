// ExactMethodForDGP.cpp: implementation of the CExactDGPMethod class.
//
//////////////////////////////////////////////////////////////////////
#include "ExactDGPMethod.h"
#include <windows.h>
#include <fstream>
#include <iterator>
#include <cassert>
#include "..\\Model3D\\Parameters.h"
using namespace std;

namespace Geodesic
{
	using namespace Model3D;
	void CExactDGPMethod::Initialize()
	{
		CDistanceApproach::Initialize();
		m_InfoAtVertices.resize(model.GetNumOfVerts());
	}

	CExactDGPMethod::CExactDGPMethod(const CRichModel& inputModel, int source) : CDistanceApproach(inputModel, source)
	{
		m_nameOfAlgorithm = "Exact";
	}

	CExactDGPMethod::CExactDGPMethod(const CRichModel& inputModel, const set<int> &indexOfSourceVerts) : CDistanceApproach(inputModel, indexOfSourceVerts)
	{
		m_nameOfAlgorithm = "Exact";
	}

	CExactDGPMethod::CExactDGPMethod(const CRichModel& inputModel, const set<int> &indexOfSourceVerts, const set<int> &destinations) : CDistanceApproach(inputModel, indexOfSourceVerts, destinations)
	{
		m_nameOfAlgorithm = "Exact";
	}

	CExactDGPMethod::CExactDGPMethod(const CRichModel& inputModel, const map<int, double> &indexOfSourceVerts) : CDistanceApproach(inputModel, indexOfSourceVerts)
	{
		m_nameOfAlgorithm = "Exact";
	}

	CExactDGPMethod::CExactDGPMethod(const CRichModel& inputModel, const map<int, double> &indexOfSourceVerts, const set<int> &destinations) : CDistanceApproach(inputModel, indexOfSourceVerts, destinations)
	{
		m_nameOfAlgorithm = "Exact";
	}

	vector<EdgePoint> CExactDGPMethod::BacktraceShortestPath(int end) const
	{
		if (m_InfoAtVertices[end].birthTimeForCheckingValidity == -1
			|| m_InfoAtVertices[end].indexOfDirectParent == -1)
		{
			assert(model.GetNumOfComponents() != 1 || model.Neigh(end).empty());
			return vector<EdgePoint>();
		}
		vector<EdgePoint> path;
		vector<int> vertexNodes;
		int index = end;
		vertexNodes.push_back(index);
		while (m_InfoAtVertices[index].indexOfDirectParent != -1)
		{
			int indexOfParent = m_InfoAtVertices[index].indexOfDirectParent;
			if (m_InfoAtVertices[index].fParentIsPseudoSource)
			{
				index = indexOfParent;
			}
			else
			{
				index = m_InfoAtVertices[index].indexOfRootVertOfDirectParent;
			}
			vertexNodes.push_back(index);
		};
		int indexOfSourceVert = index;

		for (int i = 0; i < (int)vertexNodes.size() - 1; ++i)
		{
			int lastVert = vertexNodes[i];
			path.push_back(EdgePoint(lastVert));
			if (m_InfoAtVertices[lastVert].fParentIsPseudoSource)
			{
				continue;
			}
			int parentEdgeIndex = m_InfoAtVertices[lastVert].indexOfDirectParent;
			int edgeIndex = model.Edge(parentEdgeIndex).indexOfReverseEdge;
			Eigen::Vector2d coord(model.GetNew2DCoordinatesByReversingCurrentEdge(parentEdgeIndex, model.Edge(parentEdgeIndex).coordOfOppositeVert));

			double proportion = 1 - m_InfoAtVertices[lastVert].entryProp;
			while (true)
			{
				path.push_back(EdgePoint(edgeIndex, proportion));
				if (model.Edge(edgeIndex).indexOfOppositeVert == vertexNodes[i + 1])
					break;
				double oldProprotion = proportion;
				proportion = model.ProportionOnLeftEdgeByImage(edgeIndex, coord, oldProprotion);
				if (abs(proportion - 1) < 1e-2)
				{
					vector<EdgePoint> path2 = BacktraceShortestPath(model.Edge(edgeIndex).indexOfOppositeVert);
					reverse(path.begin(), path.end());
					copy(path.begin(), path.end(), back_inserter(path2));
					return path2;
				}
				else if (proportion >= 0 && proportion <= 1)
				{
					proportion = max(proportion, 0);
					coord = model.GetNew2DCoordinatesByRotatingAroundLeftChildEdge(edgeIndex, coord);
					edgeIndex = model.Edge(edgeIndex).indexOfLeftEdge;
					//rightLen = disToAngle;				
				}
				else
				{
					proportion = model.ProportionOnRightEdgeByImage(edgeIndex, coord, oldProprotion);
					proportion = max(proportion, 0);
					proportion = min(proportion, 1);
					coord = model.GetNew2DCoordinatesByRotatingAroundRightChildEdge(edgeIndex, coord);
					edgeIndex = model.Edge(edgeIndex).indexOfRightEdge;
				}
			};
		}
		path.push_back(EdgePoint(indexOfSourceVert));
		reverse(path.begin(), path.end());
		return path;
	}

	int CExactDGPMethod::GetAncestor(int vertex) const
	{
		return m_InfoAtVertices[vertex].indexOfAncestor;
	}

	void CExactDGPMethod::CollectExperimentalResults()
	{
		m_memory = ((double)model.GetNumOfVerts() * sizeof InfoAtVertex) / 1024 / 1024;
		for (int i = 0; i < m_scalarField.size(); ++i)
		{
			m_scalarField[i] = m_InfoAtVertices[i].disUptodate;
		}
		CDistanceApproach::CollectExperimentalResults();
	}

	void CExactDGPMethod::Dispose()
	{
		//Do nothing...
	}
}