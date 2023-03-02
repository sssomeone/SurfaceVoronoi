#pragma once
#include "BaseModel.h"
#include <CGAL\Polygon_2.h>

namespace Model3D
{
	template<typename K>
	void TransformPolygon2Model3D(const CGAL::Polygon_2<K>& poly, CBaseModel& model)
	{
		vector<CPoint3D> vertList;
		vector<CBaseModel::CFace> faceList;
		int n = poly.size();
		for (int i = 0; i < n; ++i)
			vertList.push_back(CPoint3D(poly.vertex(i).x(), poly.vertex(i).y(), 0));
		for (int i = 0; i < n; ++i)
		{
			double len1 = sqrt((poly.vertex(i) - poly.vertex((i + n - 1) % n)).squared_length());
			double len2 = sqrt((poly.vertex(i) - poly.vertex((i + 1) % n)).squared_length());
			vertList.push_back(CPoint3D(poly.vertex(i).x(), poly.vertex(i).y(), (len1 + len2) / 2));
		}
		for (int i = 0; i < n; ++i)
		{
			faceList.push_back(CBaseModel::CFace(i, (i + 1) % n, (i + 1) % n + n));
			faceList.push_back(CBaseModel::CFace(i, (i + 1) % n + n, i + n));
		}
		model = CBaseModel(vertList, faceList);
	}

	template<typename K>
	CBaseModel TransformPolygon2Model3D(const CGAL::Polygon_2<K>& poly)
	{
		vector<CPoint3D> vertList;
		vector<CBaseModel::CFace> faceList;
		int n = poly.size();
		for (int i = 0; i < n; ++i)
			vertList.push_back(CPoint3D(poly.vertex(i).x(), poly.vertex(i).y(), 0));
		for (int i = 0; i < n; ++i)
		{
			double len1 = sqrt((poly.vertex(i) - poly.vertex((i + n - 1) % n)).squared_length());
			double len2 = sqrt((poly.vertex(i) - poly.vertex((i + 1) % n)).squared_length());
			vertList.push_back(CPoint3D(poly.vertex(i).x(), poly.vertex(i).y(), (len1 + len2) / 2));
		}
		for (int i = 0; i < n; ++i)
		{
			faceList.push_back(CBaseModel::CFace(i, (i + 1) % n, (i + 1) % n + n));
			faceList.push_back(CBaseModel::CFace(i, (i + 1) % n + n, i + n));
		}
		return CBaseModel(vertList, faceList);
	}
}