#include "IntrinsicSDF.h"
#include <algorithm>
#include <iterator>
using namespace std;
namespace Geodesic
{
	using namespace Model3D;
	CIntrinsicSDF::CIntrinsicSDF(const CRichModel& model, int source/*, const vector<double> &moduleOfEdges, const vector<double>& moduleOfVertices*/) : CXin_Wang(model, source)/*, moduleOfEdges(moduleOfEdges), moduleOfVertices(moduleOfVertices)*/
	{
		m_nameOfAlgorithm = "IntrinsicSDF";
		m_refISDF = FLT_MAX;
	}

	CIntrinsicSDF::CIntrinsicSDF(const CRichModel& model, int source, double refISDF) : CXin_Wang(model, source)/*, moduleOfEdges(moduleOfEdges), moduleOfVertices(moduleOfVertices)*/
	{
		m_nameOfAlgorithm = "IntrinsicSDF";
		m_refISDF = refISDF;
	}

	void CIntrinsicSDF::ComputeModuleOfVerticesAndEdges(const CRichModel& model, vector<double>& moduleOfEdges, vector<double>& moduleOfVertices)
	{
		moduleOfVertices.resize(model.GetNumOfVerts(), 0);
		moduleOfEdges.resize(model.GetNumOfEdges(), 0);

		for (int i = 0; i < model.GetNumOfVerts(); ++i)
		{
			double maxLen = 0;
			for (int j = 0; j < model.Neigh(i).size(); ++j)
			{
				if (maxLen < model.Edge(model.Neigh(i)[j].first).length)
					maxLen = model.Edge(model.Neigh(i)[j].first).length;
			}
			moduleOfVertices[i] = maxLen;
		}
		for (int i = 0; i < model.GetNumOfEdges(); ++i)
		{
			if (model.IsExtremeEdge(i))
				continue;
			double maxLen = 0;
			if (maxLen < model.Edge(model.Edge(i).indexOfLeftEdge).length)
				maxLen = model.Edge(model.Edge(i).indexOfLeftEdge).length;
			if (maxLen < model.Edge(model.Edge(i).indexOfRightEdge).length)
				maxLen = model.Edge(model.Edge(i).indexOfRightEdge).length;
			int reverseEdge = model.Edge(i).indexOfReverseEdge;
			if (model.IsExtremeEdge(reverseEdge))
				continue;
			if (maxLen < model.Edge(model.Edge(reverseEdge).indexOfLeftEdge).length)
				maxLen = model.Edge(model.Edge(reverseEdge).indexOfLeftEdge).length;
			if (maxLen < model.Edge(model.Edge(reverseEdge).indexOfRightEdge).length)
				maxLen = model.Edge(model.Edge(reverseEdge).indexOfRightEdge).length;
			moduleOfEdges[i] = maxLen;
		}
	}

	void CIntrinsicSDF::Initialize()
	{
		CChen_Han::Initialize();
		m_windowsOnEdges.resize(model.GetNumOfEdges());
		m_wavefront = 0;
		m_maxModuleOfSweptEdges = 0;
		m_nWindowsKeptOnEdges = 0;
		//m_numOfSweptSaddleVertices = 0;
		m_numOfSweptVertices = 0;
		m_nTimesForComputingLoops = 0;
		m_bestWindowPair.loopLength = DBL_MAX;
	}

	void CIntrinsicSDF::Dispose()
	{
		m_windowsOnEdges.clear();
		//	m_windowsAroundSaddleVertices.clear();
		CXin_Wang::Dispose();
	}

	void CIntrinsicSDF::Propagate()
	{
		int terminationTimes(0);
		ComputeChildrenOfSource();
		bool fFromQueueOfPseudoSources = UpdateTreeDepthBackWithChoice();
		while (!m_QueueForPseudoSources.empty() || !m_QueueForWindows.empty())
		{
			if (++terminationTimes > 1000 * model.GetNumOfVerts())
				return;
			if (m_QueueForWindows.size() > m_nMaxLenOfWindowQueue)
				m_nMaxLenOfWindowQueue = (int)m_QueueForWindows.size();
			if (m_QueueForPseudoSources.size() > m_nMaxLenOfPseudoSourceQueue)
				m_nMaxLenOfPseudoSourceQueue = m_QueueForPseudoSources.size();
			if (m_QueueForWindows.size() + m_QueueForPseudoSources.size() > m_maxLenOfQueue)
				m_maxLenOfQueue = m_QueueForWindows.size() + m_QueueForPseudoSources.size();
			if (fFromQueueOfPseudoSources) //pseudosource
			{
				int indexOfVert = m_QueueForPseudoSources.top().indexOfVert;
				{
					//if (m_maxModuleOfSweptEdges < moduleOfVertices[indexOfVert])
					//	m_maxModuleOfSweptEdges = moduleOfVertices[indexOfVert];
					m_wavefront = m_QueueForPseudoSources.top().disUptodate;

					FigureOutNewLoop(indexOfVert);
					m_numOfSweptVertices++;

					if (model.isBoundaryVert(indexOfVert))
						return;
					if (ISDFisDetermined())
						return;
					//if (min(m_bestWindowPair.loopLength, 2 * m_wavefront) > m_refISDF + 1.5 * model.m_averageEdgeLength)
					//{
					//	return;
					//}
				}
				m_QueueForPseudoSources.pop();
				ComputeChildrenOfPseudoSource(indexOfVert);
			}
			else
			{
				QuoteWindow quoteW = m_QueueForWindows.top();
				{
					m_wavefront = m_QueueForWindows.top().disUptodate;
					if (model.IsExtremeEdge(quoteW.pWindow->indexOfCurEdge))
						return;
					//all windows come here...
					int edge1 = model.Edge(quoteW.pWindow->indexOfCurEdge).indexOfReverseEdge;
					int edge2 = model.Edge(edge1).indexOfLeftEdge;
					int edge3 = model.Edge(edge1).indexOfRightEdge;
					if (m_maxModuleOfSweptEdges
						< model.Edge(edge1).length)
					{
						m_maxModuleOfSweptEdges = model.Edge(edge1).length;
					}
					if (m_maxModuleOfSweptEdges
						< model.Edge(edge2).length)
					{
						m_maxModuleOfSweptEdges = model.Edge(edge2).length;
					}
					if (m_maxModuleOfSweptEdges
						< model.Edge(edge3).length)
					{
						m_maxModuleOfSweptEdges = model.Edge(edge3).length;
					}
					BalancedWindow bw(*quoteW.pWindow);
					bw.minDis = quoteW.disUptodate;
					const CRichModel::CEdge& edge = model.Edge(quoteW.pWindow->indexOfCurEdge);
					//int leftVert = edge.indexOfLeftVert;
					double detaX = quoteW.pWindow->coordOfPseudoSource(0) - quoteW.pWindow->proportions[1] * edge.length;
					double rightLen = sqrt(detaX * detaX + quoteW.pWindow->coordOfPseudoSource(1) * quoteW.pWindow->coordOfPseudoSource(1));
					//int rightVert = edge.indexOfRightVert;
					detaX = quoteW.pWindow->coordOfPseudoSource(0) - quoteW.pWindow->proportions[0] * edge.length;
					double leftLen = sqrt(detaX * detaX + quoteW.pWindow->coordOfPseudoSource(1) * quoteW.pWindow->coordOfPseudoSource(1));
					bw.maxDis = quoteW.pWindow->disToRoot + max(leftLen, rightLen);

					FigureOutNewLoop(bw);

					m_windowsOnEdges[quoteW.pWindow->indexOfCurEdge].push_back(bw);

					if (ISDFisDetermined())
						return;
					//if (min(m_bestWindowPair.loopLength, 2 * m_wavefront) > m_refISDF + 1.5 * model.m_averageEdgeLength)
					//{
					//	return;
					//}
				}
				m_QueueForWindows.pop();
				ComputeChildrenOfWindow(quoteW);
				delete quoteW.pWindow;
			}
			fFromQueueOfPseudoSources = UpdateTreeDepthBackWithChoice();
		}
	}

	void CIntrinsicSDF::CollectExperimentalResults()
	{
		m_nWindowsKeptOnEdges = 0;
		for (int i = 0; i < m_windowsOnEdges.size(); ++i)
		{
			for (list<BalancedWindow>::const_iterator it2 = m_windowsOnEdges[i].begin();
				it2 != m_windowsOnEdges[i].end(); ++it2)
			{
				if (!CannotDetermineShortestLoop(*it2))
					m_nWindowsKeptOnEdges++;
			}
		}
		m_windowsOnEdges.clear();


		m_memory = ((double)m_numOfSweptVertices * sizeof InfoAtVertex
			+ (double)m_nWindowsKeptOnEdges * sizeof BalancedWindow
			+ (double)m_maxLenOfQueue * sizeof(Window *)) / 1024 / 1024;
		for (int i = 0; i < m_scalarField.size(); ++i)
		{
			m_scalarField[i] = m_InfoAtVertices[i].disUptodate;
		}
		CDistanceApproach::CollectExperimentalResults();
	}

	vector<EdgePoint> CIntrinsicSDF::BacktraceShortestLoop() const
	{
		if (m_bestWindowPair.loopLength > FLT_MAX)
			return vector<EdgePoint>();
		//cerr << dec;
		//cerr.precision(5);
		//cerr.setf(ios::fixed);
		//cerr << "LoopLength = " << m_bestWindowPair.loopLength << endl;
		//cerr << "isVertex1 = " << m_bestWindowPair.entryPoint1.isVertex << endl;
		//cerr << "isVertex2 = " << m_bestWindowPair.entryPoint2.isVertex << endl;
		vector<EdgePoint> path1 = BacktraceShortestPath(m_bestWindowPair.meetingPoint, m_bestWindowPair.entryPoint1);

		vector<EdgePoint> path2 = BacktraceShortestPath(m_bestWindowPair.meetingPoint, m_bestWindowPair.entryPoint2);

		reverse(path2.begin(), path2.end());
		//#include <iterator>
		copy(path2.begin() + 1, path2.end(), back_inserter(path1));
		return path1;
	}

	bool CIntrinsicSDF::ISDFisDetermined() const
	{
		if (m_bestWindowPair.loopLength > FLT_MAX)
			return false;
		//please consider a special case:
		//a window from a vertex meet another window
		return m_wavefront > m_bestWindowPair.loopLength / 2 + m_maxModuleOfSweptEdges;
	}

	double CIntrinsicSDF::GetIntrinsicSDF() const
	{
		//if (m_sources.begin()->first == 1985)
		//{
		//	cerr << "loopLength 1985: " << m_bestWindowPair.loopLength << endl;
		//	cerr << "ref 1985: " << m_refISDF << endl;
		//}
		//if (min(m_bestWindowPair.loopLength, 2 * m_wavefront) > m_refISDF + 1.5 * model.m_averageEdgeLength)
		//{
		//	return 0;
		//}	
		if (m_bestWindowPair.loopLength < FLT_MAX)
			return m_bestWindowPair.loopLength;

		return 2 * m_wavefront;
	}

	//void CIntrinsicSDF::AddIntoQueueOfWindows(QuoteWindow& quoteW)
	//{
	//	quoteW.disUptodate = GetMinDisOfWindow(*quoteW.pWindow);
	//		
	//	if (!CheckValidityWithXinWangFiltering(*quoteW.pWindow))
	//	{
	//		delete quoteW.pWindow;
	//		return;
	//	}
	//	m_QueueForWindows.push(quoteW);
	//	++m_nCountOfWindows;
	//}
	//
	//void CIntrinsicSDF::AddIntoQueueOfPseudoSources(QuoteInfoAtVertex& quoteOfPseudoSource)
	//{	
	//	m_QueueForPseudoSources.push(quoteOfPseudoSource);
	//}

	void CIntrinsicSDF::FigureOutNewLoop(const BalancedWindow& w)
	{
		//window and window meet at edge
		int reverseEdge = model.Edge(w.indexOfCurEdge).indexOfReverseEdge;
		//if (w.indexOfRootVertex == model.Edge(reverseEdge).indexOfOppositeVert
		//	|| model.Edge(reverseEdge).angleOpposite > M_PI - 6.0 * M_PI / 180)
		//{
		//	//error prone
		//	return;
		//}
		list<BalancedWindow>& wList = m_windowsOnEdges[reverseEdge];
		Eigen::Vector2d rootOfW = model.GetNew2DCoordinatesByReversingCurrentEdge(w.indexOfCurEdge, w.coordOfPseudoSource);

		for (list<BalancedWindow>::iterator it = wList.begin();
			it != wList.end(); ++it)
		{
			if (CannotDetermineShortestLoop(*it))
			{
				if (CanBeRemoved(*it))
				{
					list<BalancedWindow>::iterator oldIt = it;
					++it;
					wList.erase(oldIt);
					--it;
				}
				continue;
			}
			//important
			double leftProp = max(it->proportions[0], 1 - w.proportions[1]);
			double rightProp = min(it->proportions[1], 1 - w.proportions[0]);
			if (leftProp <= rightProp)
			{
				double detaX = rootOfW(0) - it->coordOfPseudoSource(0);
				double detaY = rootOfW(1) - it->coordOfPseudoSource(1);
				double disSum = w.disToRoot + it->disToRoot + sqrt(detaX * detaX + detaY * detaY);
				double prop = model.ProportionOnEdgeByImage(reverseEdge, rootOfW(0), rootOfW(1), it->coordOfPseudoSource(0), it->coordOfPseudoSource(1));
				if (disSum < m_bestWindowPair.loopLength
					&& prop >= leftProp - 1e-4 && prop <= rightProp + 1e-4)
				{
					m_nTimesForComputingLoops++;
					m_bestWindowPair.loopLength = disSum;
					m_bestWindowPair.meetingPoint.isVertex = false;
					m_bestWindowPair.meetingPoint.index = reverseEdge;
					m_bestWindowPair.meetingPoint.proportion
						= prop;
					if (w.fDirectParenIsPseudoSource)
					{
						m_bestWindowPair.entryPoint1.isVertex = true;
						m_bestWindowPair.entryPoint1.index = w.indexOfRootVertex;
					}
					else
					{
						m_bestWindowPair.entryPoint1.isVertex = false;
						m_bestWindowPair.entryPoint1.index = w.GetDirectParentEdge(model);

						if (model.Edge(m_bestWindowPair.entryPoint1.index).indexOfLeftEdge == w.indexOfCurEdge)
						{
							m_bestWindowPair.entryPoint1.proportion = 1 - model.ProportionOnRightEdgeByImage(reverseEdge,
								it->coordOfPseudoSource,
								m_bestWindowPair.meetingPoint.proportion);
						}
						else
						{
							assert(model.Edge(m_bestWindowPair.entryPoint1.index).indexOfRightEdge == w.indexOfCurEdge);
							m_bestWindowPair.entryPoint1.proportion = 1 - model.ProportionOnLeftEdgeByImage(reverseEdge,
								it->coordOfPseudoSource,
								m_bestWindowPair.meetingPoint.proportion);
						}
					}

					if (it->fDirectParenIsPseudoSource)
					{
						m_bestWindowPair.entryPoint2.isVertex = true;
						m_bestWindowPair.entryPoint2.index = it->indexOfRootVertex;
					}
					else
					{
						m_bestWindowPair.entryPoint2.isVertex = false;
						m_bestWindowPair.entryPoint2.index = it->GetDirectParentEdge(model);
						if (model.Edge(m_bestWindowPair.entryPoint2.index).indexOfLeftEdge == reverseEdge)
						{
							//fix bug...
							m_bestWindowPair.entryPoint2.proportion = 1 - model.ProportionOnRightEdgeByImage(w.indexOfCurEdge,
								w.coordOfPseudoSource,
								1 - m_bestWindowPair.meetingPoint.proportion);
						}
						else
						{
							assert(model.Edge(it->GetDirectParentEdge(model)).indexOfRightEdge == reverseEdge);
							//fix bug...
							m_bestWindowPair.entryPoint2.proportion = 1 - model.ProportionOnLeftEdgeByImage(w.indexOfCurEdge,
								w.coordOfPseudoSource,
								1 - m_bestWindowPair.meetingPoint.proportion);
						}
					}
					m_candidateWindowPairs.push_back(m_bestWindowPair);
				}
			}
		}
	}

	void CIntrinsicSDF::FigureOutNewLoop(int indexOfUpdatedVertex)
	{
		if (model.IsStronglyConvexVert(indexOfUpdatedVertex))
			return;
		for (int i = 0; i < model.Neigh(indexOfUpdatedVertex).size(); ++i)
		{
			int neigh = model.Edge(model.Neigh(indexOfUpdatedVertex)[i].first).indexOfRightVert;
			if (CannotDetermineShortestLoop(indexOfUpdatedVertex, neigh))
				continue;
			if (model.IsStronglyConvexVert(neigh))
				continue;
			if (m_InfoAtVertices[neigh].fParentIsPseudoSource
				&& m_InfoAtVertices[neigh].indexOfDirectParent == indexOfUpdatedVertex)
				continue;
			if (m_InfoAtVertices[indexOfUpdatedVertex].fParentIsPseudoSource
				&& m_InfoAtVertices[indexOfUpdatedVertex].indexOfDirectParent == neigh)
				continue;
			double disSum = m_InfoAtVertices[neigh].disUptodate
				+ m_InfoAtVertices[indexOfUpdatedVertex].disUptodate
				+ model.Edge(model.Neigh(indexOfUpdatedVertex)[i].first).length;
			if (disSum < m_bestWindowPair.loopLength)
			{
				EdgePoint pt1 = m_InfoAtVertices[indexOfUpdatedVertex].GetEntryPoint();
				EdgePoint pt2 = m_InfoAtVertices[neigh].GetEntryPoint();

				pair<double, double> angles1 = model.GetTwoSplitAngles(indexOfUpdatedVertex, pt1, EdgePoint(neigh));
				bool flag1 = angles1.first > M_PI - 5 * M_PI / 180 && angles1.second > M_PI - 5 * M_PI / 180;
				pair<double, double> angles2 = model.GetTwoSplitAngles(neigh, EdgePoint(indexOfUpdatedVertex), pt2);
				bool flag2 = angles2.first > M_PI - 5 * M_PI / 180 && angles2.second > M_PI - 5 * M_PI / 180;
				if (flag1 && flag2)
				{
					m_nTimesForComputingLoops++;
					m_bestWindowPair.loopLength = disSum;
					m_bestWindowPair.entryPoint1 = pt1;
					m_bestWindowPair.meetingPoint = EdgePoint(indexOfUpdatedVertex);
					m_bestWindowPair.entryPoint2 = EdgePoint(neigh);
					m_candidateWindowPairs.push_back(m_bestWindowPair);
				}
			}
		}
	}

	void CIntrinsicSDF::OutputExperimentalResults() const
	{
		cout << "Experimental results are as follows:\n";
		cout << "Algorithm: " << m_nameOfAlgorithm << endl;
		cout << "Memory = " << m_memory << " Mega-bytes.\n";
		cout << "Timing = " << m_nTotalMilliSeconds << " ms.\n";
		cout << "MaxDepth = " << m_depthOfResultingTree << " levels.\n";
		cout << "MaxDis = " << GetMaxDistance() << endl;
		cout << "MaxLenOfQue = " << m_maxLenOfQueue << " elements.\n";
		cout << "TotalWindowNum = " << m_nCountOfWindows << endl;
		cout << "NumOfWindowKept = " << m_nWindowsKeptOnEdges << endl;
		cout << "TryLoopTimes = " << m_nTimesForComputingLoops << endl;
		cout << "TryTimes per window = " << (double)m_nTimesForComputingLoops / (double)m_nCountOfWindows << endl;
		cout << "IntrinsicSDF = " << GetIntrinsicSDF() << endl;
	}

	vector<EdgePoint> CIntrinsicSDF::BacktraceShortestPath(EdgePoint start, EdgePoint pivot) const
	{
		vector<EdgePoint> path;
		if (start.isVertex && start.index == m_sources.begin()->first)
		{
			path.push_back(start);
			return path;
		}
		if (pivot.isVertex)
		{
			path = CExactDGPMethod::BacktraceShortestPath(pivot.index);
			path.push_back(start);
		}
		else
		{
			path.push_back(start);

			int parentEdgeIndex = pivot.index;
			int edgeIndex = model.Edge(parentEdgeIndex).indexOfReverseEdge;
			Eigen::Vector2d coord;
			if (start.isVertex)
				coord = model.GetNew2DCoordinatesByReversingCurrentEdge(parentEdgeIndex, model.Edge(parentEdgeIndex).coordOfOppositeVert);
			else
			{
				if (model.Edge(start.index).indexOfLeftEdge != edgeIndex
					&& model.Edge(start.index).indexOfRightEdge != edgeIndex)
				{
					start = EdgePoint(model.Edge(start.index).indexOfReverseEdge, 1 - start.proportion);
				}

				assert(model.Edge(start.index).indexOfLeftEdge == edgeIndex
					|| model.Edge(start.index).indexOfRightEdge == edgeIndex);

				coord(0) = start.proportion * model.Edge(start.index).length;
				coord(1) = 0;
				if (model.Edge(start.index).indexOfLeftEdge == edgeIndex)
					coord = model.GetNew2DCoordinatesByRotatingAroundLeftChildEdge(start.index, coord);
				else
					coord = model.GetNew2DCoordinatesByRotatingAroundRightChildEdge(start.index, coord);
			}

			double proportion = 1 - pivot.proportion;
			while (true)
			{
				path.push_back(EdgePoint(edgeIndex, proportion));

				double oldProprotion = proportion;
				proportion = model.ProportionOnLeftEdgeByImage(edgeIndex, coord, oldProprotion);
				if (model.Edge(edgeIndex).indexOfOppositeVert == m_sources.begin()->first)
				{
					path.push_back(EdgePoint(m_sources.begin()->first));
					reverse(path.begin(), path.end());
					return path;
				}
				else if (abs(proportion - 1) < 1e-2)
				{
					vector<EdgePoint> path2 = CExactDGPMethod::BacktraceShortestPath(model.Edge(edgeIndex).indexOfOppositeVert);
					reverse(path.begin(), path.end());
					copy(path.begin(), path.end(), back_inserter(path2));
					return path2;
				}
				else if (abs(proportion) < 1e-2)
				{
					vector<EdgePoint> path2 = CExactDGPMethod::BacktraceShortestPath(model.Edge(edgeIndex).indexOfLeftVert);
					path.pop_back();
					reverse(path.begin(), path.end());
					copy(path.begin(), path.end(), back_inserter(path2));
					return path2;
				}
				else if (proportion > 0 && proportion < 1)
				{
					coord = model.GetNew2DCoordinatesByRotatingAroundLeftChildEdge(edgeIndex, coord);
					edgeIndex = model.Edge(edgeIndex).indexOfLeftEdge;
				}
				else
				{
					proportion = model.ProportionOnRightEdgeByImage(edgeIndex, coord, oldProprotion);
					proportion = max(proportion, 0.);
					proportion = min(proportion, 1.);
					coord = model.GetNew2DCoordinatesByRotatingAroundRightChildEdge(edgeIndex, coord);
					edgeIndex = model.Edge(edgeIndex).indexOfRightEdge;
				}
			}
		}
		return path;
	}

	bool CIntrinsicSDF::CanBeRemoved(const BalancedWindow& w) const
	{
		int edge1 = w.indexOfCurEdge;
		double module = model.Edge(edge1).length;

		if (!model.IsExtremeEdge(edge1))
		{
			int edge2 = model.Edge(edge1).indexOfLeftEdge;
			int edge3 = model.Edge(edge1).indexOfRightEdge;

			if (module < model.Edge(edge2).length)
				module = model.Edge(edge2).length;
			if (module < model.Edge(edge3).length)
				module = model.Edge(edge3).length;
		}
		if (w.maxDis + 2 * module /*/ 2.0*/ < m_wavefront - LengthTolerance)
			return true;
		return false;
	}

	bool CIntrinsicSDF::CannotDetermineShortestLoop(const BalancedWindow& w) const
	{
		//check the windows on the back side
		//has been considered already.
		//考虑到一边可能有顶点，我们放松下面的要求
		//if (w.maxDis + 2 * moduleOfEdges[w.indexOfCurEdge] /*/ 2.0*/ < m_wavefront - LengthTolerance
		//	|| w.minDis > m_wavefront + LengthTolerance)
		//	return true;

		int edge1 = w.indexOfCurEdge;
		double module = model.Edge(edge1).length;

		if (!model.IsExtremeEdge(edge1))
		{
			int edge2 = model.Edge(edge1).indexOfLeftEdge;
			int edge3 = model.Edge(edge1).indexOfRightEdge;

			if (module < model.Edge(edge2).length)
				module = model.Edge(edge2).length;
			if (module < model.Edge(edge3).length)
				module = model.Edge(edge3).length;
		}
		if (w.maxDis + 2 * module /*/ 2.0*/ < m_wavefront - LengthTolerance
			|| w.minDis > m_wavefront + LengthTolerance)
			return true;

		//too far from the source	
		if (w.minDis > m_bestWindowPair.loopLength / 2.0 + LengthTolerance)
			return true;

		return false;
	}

	bool CIntrinsicSDF::CannotDetermineShortestLoop(int current, int before) const
	{
		//the current vertx comes later.
		if (m_InfoAtVertices[before].disUptodate > m_InfoAtVertices[current].disUptodate + LengthTolerance)
			return true;
		//the loop is longer than before
		if (m_InfoAtVertices[before].disUptodate
			+ m_InfoAtVertices[current].disUptodate
			+ model.Edge(model.GetEdgeIndexFromTwoVertices(before, current)).length
		> m_bestWindowPair.loopLength)
			return true;
		return false;
	}

	void CIntrinsicSDF::SaveAllLoops(const string& filepath) const
	{
		for (int i = 0; i < m_candidateWindowPairs.size(); ++i)
		{
			vector<EdgePoint> path = BacktraceShortestPath(m_candidateWindowPairs[i].meetingPoint, m_candidateWindowPairs[i].entryPoint1);
			vector<EdgePoint> path2 = BacktraceShortestPath(m_candidateWindowPairs[i].meetingPoint, m_candidateWindowPairs[i].entryPoint2);

			reverse(path2.begin(), path2.end());
			//#include <iterator>
			copy(path2.begin() + 1, path2.end(), back_inserter(path));
			char filename[256];
			sprintf_s(filename, "%s_%04d.obj", filepath.c_str(), i + 1);
			model.SavePathToObj(path, filename);
		}
	}

	double CIntrinsicSDF::DistanceBetweenTwoLoops(const vector<EdgePoint>& loop1, const vector<EdgePoint>& loop2, const CRichModel& model)
	{
		double finalDis = 0;
		for (int i = 0; i < loop1.size(); ++i)
		{
			double dis = FLT_MAX;
			for (int j = 0; j < loop2.size(); ++j)
			{
				double tmpDis = (loop1[i].Get3DPoint(model) - loop2[j].Get3DPoint(model)).Len();
				if (tmpDis < dis)
					dis = tmpDis;
			}
			if (dis > finalDis)
				finalDis = dis;
		}
		for (int i = 0; i < loop2.size(); ++i)
		{
			double dis = FLT_MAX;
			for (int j = 0; j < loop1.size(); ++j)
			{
				double tmpDis = (loop1[j].Get3DPoint(model) - loop2[i].Get3DPoint(model)).Len();
				if (tmpDis < dis)
					dis = tmpDis;
			}
			if (dis > finalDis)
				finalDis = dis;
		}
		return finalDis;
	}

	bool CIntrinsicSDF::IsLocatedOnCorner() const
	{
		if (min(m_bestWindowPair.loopLength, 2 * m_wavefront) > m_refISDF + 1.5 * model.m_maxEdgeLength)
		{
			return true;
		}
		return false;
	}
}