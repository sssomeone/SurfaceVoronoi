#pragma once
#include "..\\Model3D\\EdgePoint.h"
#include "..\\Model3D\\RichModel.h"

//#if defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Debug\\Model3D.lib")
//#endif
//#if !defined(_DEBUG) && defined(_WIN64)
//#pragma comment(lib, "..\\x64\\Release\\Model3D.lib")
//#endif
namespace Geodesic
{
	using namespace Model3D;
	class CDistanceApproach
	{
	protected:
		vector<double> m_scalarField;
		double m_maxDisValue;
		const CRichModel& model;
		map<int, double> m_sources;
		set<int> m_destinations;
		string m_nameOfAlgorithm;
		long m_maxLenOfQueue;
		long m_depthOfResultingTree;

		long m_nTotalMilliSeconds;
	protected:
		virtual void Initialize();
		virtual void Dispose();
		virtual void Propagate() = 0;
		virtual void CollectExperimentalResults();
	public:
		double m_memory;
		CDistanceApproach(const CRichModel& model);
		CDistanceApproach(const CRichModel& model, int source);
		CDistanceApproach(const CRichModel& model, const map<int, double>& sources);
		CDistanceApproach(const CRichModel& model, const map<int, double>& sources, const set<int> &destinations);
		CDistanceApproach(const CRichModel& model, const set<int>& sources);
		CDistanceApproach(const CRichModel& model, const set<int>& sources, const set<int>& destinations);
	public:
		virtual void Execute();
		virtual vector<EdgePoint> BacktraceShortestPath(int end) const = 0;
		vector<EdgePoint> BacktraceIsoline(double val) const;
		virtual int GetAncestor(int vIndex) const = 0;
		double GetMaxDistance() const;
		virtual  __int64 GetMaxLenOfQueue() const { return m_maxLenOfQueue; }
		virtual  __int64 GetMaxPropagationLevels()const { return m_depthOfResultingTree; }
		string GetAlgorithmName() const;
		__int64 GetRunTime() const { return m_nTotalMilliSeconds; }
		double GetMemoryCost() const { return m_memory; }
		virtual void OutputExperimentalResults() const;
		const vector<double>& GetDistanceField() const;
		//vector<double> GetNormalizedDistanceField() const;
		//vector<double> GetDistanceFieldDividedBy(double denominator) const;
		static vector<double> DiffDistanceField(const vector<double>& field1,
			const vector<double>& field2);
	};

	template<class T>
	vector<EdgePoint> GetShortestPathBetween(const CRichModel& model, int source, int end)
	{
		set<int> sources;
		sources.insert(source);
		set<int> destinations;
		destinations.insert(end);

		T alg(model, sources, destinations);
		alg.Execute();
		return alg.BacktraceShortestPath(end);
	}
}