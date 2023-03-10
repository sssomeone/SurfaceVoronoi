// ImprovedCHWithEdgeValve.h: interface for the CICH_WindowFiltering class.
//
//////////////////////////////////////////////////////////////////////
#pragma once

#include "Chen_Han.h"
namespace Geodesic
{
	using namespace Model3D;
	class CICH_WindowFiltering : public CChen_Han
	{
	protected:
		virtual bool CheckValidityWithXinWangFiltering(Window& w) const;
		virtual void AddIntoQueueOfWindows(QuoteWindow& quoteW);
	public:
		CICH_WindowFiltering(const CRichModel& model, int source);
		CICH_WindowFiltering(const CRichModel& model, const map<int, double>& sources);
		CICH_WindowFiltering(const CRichModel& model, const map<int, double>& sources, const set<int> &destinations);
		CICH_WindowFiltering(const CRichModel& model, const set<int>& sources);
		CICH_WindowFiltering(const CRichModel& model, const set<int>& sources, const set<int>& destinations);
	};
}