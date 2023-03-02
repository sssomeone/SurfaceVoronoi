#pragma once
#include "xin_wang.h"
namespace Geodesic
{
	using namespace Model3D;
	class CIntrinsicSDF : public CXin_Wang
	{
	protected:
		//const vector<double>& moduleOfEdges;
		//const vector<double>& moduleOfVertices;
		double m_refISDF;
		struct BalancedWindow : public Window
		{
			double minDis;
			double maxDis;
			BalancedWindow() {}
			BalancedWindow(const Window& w)
			{
				this->fIsOnLeftSubtree = w.fIsOnLeftSubtree;
				this->fBrachParentIsPseudoSource = w.fBrachParentIsPseudoSource;
				this->fDirectParentEdgeOnLeft = w.fDirectParentEdgeOnLeft;
				this->fDirectParenIsPseudoSource = w.fDirectParenIsPseudoSource;
				this->birthTimeOfParent = w.birthTimeOfParent;
				this->indexOfBrachParent = w.indexOfBrachParent;
				this->indexOfRootVertex = w.indexOfRootVertex;
				this->indexOfCurEdge = w.indexOfCurEdge;
				this->levelOnSequenceTree = w.levelOnSequenceTree;
				this->indexOfAncestor = w.indexOfAncestor;
				this->disToRoot = w.disToRoot;
				this->proportions[0] = w.proportions[0];
				this->proportions[1] = w.proportions[1];
				this->entryPropOfParent = w.entryPropOfParent;
				this->coordOfPseudoSource = w.coordOfPseudoSource;
			}

			bool operator<(const BalancedWindow& other) const
			{
				if (maxDis < other.maxDis)
					return true;
				return false;
			}
		};

		struct WindowPair
		{
			EdgePoint meetingPoint;
			EdgePoint entryPoint1, entryPoint2;
			double loopLength;
			WindowPair()
			{
				loopLength = DBL_MAX;
			}
			bool operator<(const WindowPair& other) const
			{
				return loopLength < other.loopLength;
			}
		};

		vector<list<BalancedWindow>> m_windowsOnEdges;

		double m_wavefront;//done
		double m_maxModuleOfSweptEdges;//done

		__int64 m_numOfSweptVertices;//done
		__int64 m_nTimesForComputingLoops;
	protected:
		bool ISDFisDetermined() const;//done
		bool CannotDetermineShortestLoop(const BalancedWindow& w) const;
		bool CanBeRemoved(const BalancedWindow& w) const;
		bool CannotDetermineShortestLoop(int current, int before) const;
	protected:
		virtual void Initialize();//done
		virtual void Dispose();//done
		virtual void Propagate();
		virtual void CollectExperimentalResults();//done
		virtual void FigureOutNewLoop(const BalancedWindow& w);//done
		virtual void FigureOutNewLoop(int indexOfUpdatedVertex);//done
	public:
		__int64 m_nWindowsKeptOnEdges;//done
		WindowPair m_bestWindowPair;
		bool IsLocatedOnCorner() const;
		vector<WindowPair> m_candidateWindowPairs;
		//done
		CIntrinsicSDF(const CRichModel& model, int source);
		CIntrinsicSDF(const CRichModel& model, int source, double refISDF);
		//done
		static void ComputeModuleOfVerticesAndEdges(const CRichModel& model, vector<double>& moduleOfEdges, vector<double>& moduleOfVertices);
		vector<EdgePoint> BacktraceShortestLoop() const;//done
		void SaveAllLoops(const string& path) const;
		vector<EdgePoint> BacktraceShortestPath(EdgePoint start, EdgePoint pivot) const;//done
		double GetIntrinsicSDF() const;//done
		virtual void OutputExperimentalResults() const;//done
		static double DistanceBetweenTwoLoops(const vector<EdgePoint>& loop1, const vector<EdgePoint>& loop2, const CRichModel& model);
	};

}