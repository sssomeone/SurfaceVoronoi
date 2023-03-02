#pragma comment(lib, "libmpfr-4.lib")
#pragma comment(lib, "libgmp-10.lib")

#include "RichModel.h"
#include <iostream>
#include "CDT.h"

#include "Polygon2Model3D.hpp"

using namespace std;
using namespace Model3D;
int main()
{
	CRichModel model("..\\data\\sphere.obj");
	model.LoadModel();
	model.PrintInfo(cout);

	CGAL::Polygon_2<K> poly;
	poly.push_back(CGAL::Point_2<K>(0, 0));
	poly.push_back(CGAL::Point_2<K>(1, 0));
	poly.push_back(CGAL::Point_2<K>(1, 1));
	poly.push_back(CGAL::Point_2<K>(0, 1));
	TransformPolygon2Model3D<K>(poly, model);
	model.SaveObjFile("test.obj");

    return 0;
}

