#pragma once
#include "RichModel.h"
#include <fstream>
#include <Eigen\dense>
using namespace std;
struct FacePoint
{
	int faceID;
	double lamda1, lamda2, lamda3;

	FacePoint(){}
	FacePoint(int faceID, double lamda1, double lamda2)
	: faceID(faceID), lamda1(lamda1), lamda2(lamda2) 
	{
		lamda3 = 1 - lamda1 - lamda2;
	}
	FacePoint(const CRichModel& model, int faceID, const CPoint3D& pt)
		: faceID(faceID)
	{
		Eigen::MatrixXd M(4, 3);
		M << 1, 1, 1,
			model.Vert(model.Face(faceID)[0]).x, model.Vert(model.Face(faceID)[1]).x, model.Vert(model.Face(faceID)[2]).x,
			model.Vert(model.Face(faceID)[0]).y, model.Vert(model.Face(faceID)[1]).y, model.Vert(model.Face(faceID)[2]).y,
			model.Vert(model.Face(faceID)[0]).z, model.Vert(model.Face(faceID)[1]).z, model.Vert(model.Face(faceID)[2]).z;
		Eigen::VectorXd b(4);
		b << 1, pt.x, pt.y, pt.z;
		b = M.transpose() * b;
		M = M.transpose() * M;
		Eigen::VectorXd x = M.inverse() * b;
		lamda1 = x(0);
		lamda2 = x(1);
		lamda3 = x(2);
	}

	CPoint3D Get3DPoint(const CRichModel& model) const
	{
		return lamda1 * model.Vert(model.Face(faceID)[0])
			+ lamda2 * model.Vert(model.Face(faceID)[1])
			+ lamda3 * model.Vert(model.Face(faceID)[2]);
	}

	CPoint3D Get3DPoint_Shift(const CRichModel& model) const
	{
		return lamda1 * model.GetShiftVertex(model.Face(faceID)[0])
			+ lamda2 * model.GetShiftVertex(model.Face(faceID)[1])
			+ lamda3 * model.GetShiftVertex(model.Face(faceID)[2]);
	}

	pair<double, double> Get2DCoord(const CRichModel& model, int edgeIndex) const
	{
		double x, y;
		if (model.Edge(edgeIndex).indexOfLeftVert == model.Face(faceID)[0])
		{
			x = lamda1 * 0 + lamda2 * model.Edge(edgeIndex).length + lamda3 * model.Edge(edgeIndex).coordOfOppositeVert.first;
			y = lamda1 * 0 + lamda2 * 0 + lamda3 * model.Edge(edgeIndex).coordOfOppositeVert.second;
		}
		else if (model.Edge(edgeIndex).indexOfLeftVert == model.Face(faceID)[1])
		{
			x = lamda2 * 0 + lamda3 * model.Edge(edgeIndex).length + lamda1 * model.Edge(edgeIndex).coordOfOppositeVert.first;
			y = lamda2 * 0 + lamda3 * 0 + lamda1 * model.Edge(edgeIndex).coordOfOppositeVert.second;
		}
		else
		{
			x = lamda3 * 0 + lamda1 * model.Edge(edgeIndex).length + lamda2 * model.Edge(edgeIndex).coordOfOppositeVert.first;
			y = lamda3 * 0 + lamda1 * 0 + lamda2 * model.Edge(edgeIndex).coordOfOppositeVert.second;
		}
		return make_pair(x, y);
	}
};