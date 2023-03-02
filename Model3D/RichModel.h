// RichModel.h: interface for the CRichModel class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "BaseModel.h"
#include <cassert>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Eigen\core>
using namespace std;
namespace Model3D
{
	struct EdgePoint;
	class CRichModel : public CBaseModel
	{
	public:
		struct CEdge
		{
			int indexOfLeftVert;
			int indexOfRightVert;
			int indexOfOppositeVert;
			int indexOfLeftEdge;
			int indexOfRightEdge;
			int indexOfReverseEdge;
			int indexOfFrontFace;
			double length;
			double angleOpposite;
			Eigen::Vector2d coordOfOppositeVert;
			// |unitX   -unitY|
			// |unitY    unitX|
			Eigen::Vector2d matrixRotatedToLeftEdge;
			Eigen::Vector2d matrixRotatedToRightEdge;
			CEdge()
			{
				indexOfOppositeVert = -1;	//key	
				indexOfLeftEdge = -1;
				indexOfRightEdge = -1;
				indexOfFrontFace = -1;
				angleOpposite = 2 * M_PI;
			}
		};

	protected:
		void CreateEdgesFromVertsAndFaces();
		void CollectAndArrangeNeighs();
		void ComputeAnglesAroundVerts();
		void ComputePlanarCoordsOfIncidentVertForEdges();
		void ComputeNumOfHoles();
		void ComputeNumOfComponents();
		void ExtractBoundary();
	public:
		CRichModel(const CBaseModel &model);
		CRichModel(const string &filename);
		CRichModel(const vector<CPoint3D> &verts, const vector<CBaseModel::CFace> &faces);
		void LoadModel();
		map<CPoint3D, int> AddFaceInteriorPoints(const vector<pair<int, CPoint3D>> &pts);
		int AddFaceInteriorPoints(int faceIndex, CPoint3D pt);
		void RemoveUnreferencedVertices();
		void MakeOpen2Closed();
		void PreprocessBaseModelIntoRichModel();
		int SplitEdge(const EdgePoint& ep);
		void SplitEdgeSet(const set<int> & edgeSet, double tolerance);
		void SplitBasedOnScalarField_into_Two_Models(const vector<double>& scalarField,
			double val,
			const string& fileWithLargerScalars,
			const string& fileWithSmallerScalars) const;
		pair<CBaseModel, CBaseModel> SplitBasedOnScalarField_into_Two_Models(const vector<double>& scalarField,
			double val = 0) const;
		void SplitBasedOnScalarField_Update(const vector<double>& scalarField,
			double val = 0);
		vector<vector<EdgePoint>> GetIsoLine(const vector<double>& scalarField,
			double val) const;

		void SavePathToObj(const vector<EdgePoint>& pl, const string& filename) const;		
		void SaveIsolineToObj(const vector<EdgePoint>& isoline, const string& filename) const;
		void SaveIsolineToObj(const vector<vector<EdgePoint>>& isoline, const string& filename) const;
		void SetEdgeLength(int leftVert, int rightVert, double newLength);
		void FinishChangingEdgeLengths();
		void PrintInfo(ostream& out) const;
		double AngleSum(int vertIndex) const;
		int GetSubindexToVert(int root, int neigh) const;
		const CEdge& Edge(int edgeIndex) const;
		const vector<pair<int, double> >& Neigh(int root) const;
		//inline double Curvature(int vertIndex) const;
		//compute the proportion by two points
		double ProportionOnEdgeByImage(int edgeIndex, const Eigen::Vector2d &coord) const;
		//compute the proportion on the left edge
		double ProportionOnLeftEdgeByImage(int edgeIndex, const Eigen::Vector2d &coord, double proportion) const;
		//compute the proportion on the right edge
		double ProportionOnRightEdgeByImage(int edgeIndex, const Eigen::Vector2d &coord, double proportion) const;
		//compute the proportion by two points
		double ProportionOnEdgeByImage(int edgeIndex, double x1, double y1, double x2, double y2) const;
		Eigen::Vector2d GetNew2DCoordinatesByRotatingAroundLeftChildEdge(int edgeIndex, const Eigen::Vector2d& input2DCoordinates) const;
		Eigen::Vector2d GetNew2DCoordinatesByRotatingAroundRightChildEdge(int edgeIndex, const Eigen::Vector2d& input2DCoordinates) const;
		Eigen::Vector2d GetNew2DCoordinatesByReversingCurrentEdge(int edgeIndex, const Eigen::Vector2d& input2DCoordinates) const;
		double DistanceToOppositeAngle(int edgeIndex, const Eigen::Vector2d& coord) const;
		double DistanceToLeftVert(int edgeIndex, const Eigen::Vector2d& coord) const;
		double DistanceToRightVert(int edgeIndex, const Eigen::Vector2d& coord) const;
		int GetNumOfEdges() const;
		int GetNumOfValidDirectedEdges() const;
		int GetNumOfTotalUndirectedEdges() const;
		int GetNumOfGenera() const;
		int GetNumOfIsolated() const;
		int GetNumOfComponents() const;
		int GetNumOfBoundries() const;
		bool IsStronglyConvexVert(int index) const;
		bool IsWeaklyConvexVert(int index) const;
		bool isBoundaryVert(int index) const;
		bool IsClosedModel() const;
		bool IsExtremeEdge(int edgeIndex) const;
		bool IsStartEdge(int edgeIndex) const;
		int GetEdgeIndexFromTwoVertices(int leftVert, int rightVert) const;
		pair<double, double> GetTwoSplitAngles(int root, EdgePoint pt1, EdgePoint pt2) const;
		int IntersectQuery(int faceID, const pair<EdgePoint, EdgePoint>& seg1, const pair<EdgePoint, EdgePoint>& seg2, EdgePoint& intersection) const;
		double GetMaxEdgeLength() const;
		vector<double> ReadScalarField(const char* file) const;
		Eigen::Vector2d Get2DCoord(int face, CPoint3D pt) const;
		CPoint3D Get3DCoord(int face, Eigen::Vector2d pt) const;
	protected:
		int m_nBoundries;
		int m_nIsolatedVerts;
		int m_nComponents;

		vector<CEdge> m_Edges;
		set<int> m_UselessEdges;
		vector<vector<pair<int, double> > > m_NeighsAndAngles;
		//<strongconvex, weakconvex>
		vector<pair<bool, bool>> m_FlagsForCheckingConvexVerts;
	public:
		vector<vector<int>> m_BoundaryCurves;
		map<int, int> m_OnWhichBoundary;
		double m_maxEdgeLength;
		vector<tuple<CPoint3D,CPoint3D, CPoint3D, CPoint3D>> m_coord_systems_faces;
	};
}