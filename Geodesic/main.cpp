#pragma comment(lib, "libmpfr-4.lib")
#pragma comment(lib, "libgmp-10.lib")

#include "Xin_Wang.h"
#include <iostream>
using namespace std;

#if defined(_DEBUG) && defined(_WIN64)
#pragma comment(lib, "..\\lib\\Debug\\Model3D.lib")
#endif
#if !defined(_DEBUG) && defined(_WIN64)
#pragma comment(lib, "..\\lib\\Release\\Model3D.lib")
#endif

using namespace Model3D;
using namespace Geodesic;
void main()
{
	CRichModel model("..\\data\\sphere.obj");
	model.LoadModel();
	model.PrintInfo(cout);

	set<int> destinations;
	CXin_Wang alg(model, 0, 0);
	alg.Execute();
	cout << alg.GetMaxDistance() << endl;
	alg.GetDistanceField();
}