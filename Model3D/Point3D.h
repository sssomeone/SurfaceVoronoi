// Point3D.h: interface for the CPoint3D class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <math.h>
#define _USE_MATH_DEFINES
#include <ostream>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;
namespace Model3D
{
	struct CPoint3D
	{
	public:
		double x, y, z;
		CPoint3D();
		CPoint3D(double x, double y, double z);
		inline CPoint3D& operator +=(const CPoint3D& pt);
		inline CPoint3D& operator -=(const CPoint3D& pt);
		inline CPoint3D& operator *=(double times);
		inline CPoint3D& operator /=(double times);
		inline CPoint3D Rotate() const;
		inline double Len() const;
		inline double LenSquare() const;
		inline void Normalize();
		bool operator==(const CPoint3D& other) const
		{
			return x == other.x && y == other.y && z == other.z;
		}
		bool operator!=(const CPoint3D& other) const
		{
			return !(x == other.x && y == other.y && z == other.z);
		}
		bool operator<(const CPoint3D& other) const
		{
			if (x < other.x)
				return true;
			else if (x > other.x)
				return false;
			else if (y < other.y)
				return true;
			else if (y > other.y)
				return false;
			else if (z < other.z)
				return true;
			else if (z > other.z)
				return false;
			return false;
		}
		bool operator>(const CPoint3D& other) const
		{
			if (x > other.x)
				return true;
			else if (x < other.x)
				return false;
			else if (y > other.y)
				return true;
			else if (y < other.y)
				return false;
			else if (z > other.z)
				return true;
			else if (z < other.z)
				return false;
			return false;
		}
		friend ostream& operator<<(ostream& out, const CPoint3D& pt)
		{
			out << pt.x << " " << pt.y << " " << pt.z;
			return out;
		}
		CPoint3D GetUnitPerpendicularDir() const
		{
			vector<pair<double, int>> components;
			components.push_back(make_pair(abs(x), 1));
			components.push_back(make_pair(abs(y), 2));
			components.push_back(make_pair(abs(z), 3));
			sort(components.begin(), components.end());
			if (components[0].second == 1)
			{
				CPoint3D dir(0, z, -y);
				dir.Normalize();
				return dir;
			}
			else if (components[0].second == 2)
			{
				CPoint3D dir(-z, 0, x);
				dir.Normalize();
				return dir;
			}
			CPoint3D dir(y, -x, 0);
			dir.Normalize();
			return dir;
		}
	};
	CPoint3D CPoint3D::Rotate() const
	{
		return CPoint3D(z, x, y);
	}

	CPoint3D& CPoint3D::operator +=(const CPoint3D& pt)
	{
		x += pt.x;
		y += pt.y;
		z += pt.z;
		return *this;
	}

	CPoint3D& CPoint3D::operator -=(const CPoint3D& pt)
	{
		x -= pt.x;
		y -= pt.y;
		z -= pt.z;
		return *this;
	}

	CPoint3D& CPoint3D::operator *=(double times)
	{
		x *= times;
		y *= times;
		z *= times;
		return *this;
	}

	CPoint3D& CPoint3D::operator /=(double times)
	{
		x /= times;
		y /= times;
		z /= times;
		return *this;
	}

	double CPoint3D::Len() const
	{
		return sqrt(x * x + y * y + z * z);
	}

	double CPoint3D::LenSquare() const
	{
		return (x * x + y * y + z * z);
	}

	void CPoint3D::Normalize()
	{
		double len = Len();
		x /= len;
		y /= len;
		z /= len;
	}

	CPoint3D operator +(const CPoint3D& pt1, const CPoint3D& pt2);
	CPoint3D operator -(const CPoint3D& pt1, const CPoint3D& pt2);
	CPoint3D operator *(const CPoint3D& pt, double times);
	CPoint3D operator /(const CPoint3D& pt, double times);
	CPoint3D operator *(double times, const CPoint3D& pt);
	CPoint3D operator *(const CPoint3D& pt1, const CPoint3D& pt2);
	CPoint3D VectorCross(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3);
	double operator ^(const CPoint3D& pt1, const CPoint3D& pt2);
	double GetTriangleArea(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3);
	double AngleBetween(const CPoint3D& pt1, const CPoint3D& pt2);
	double AngleBetween(const CPoint3D& pt1, const CPoint3D& pt2, const CPoint3D& pt3);
	void VectorCross(const float* u, const float* v, float * n);
	float VectorDot(const float* u, const float* v);
	float AngleBetween(const float* u, const float* v);
	CPoint3D CombinePointAndNormalTo(const CPoint3D& pt, const CPoint3D& normal);
	CPoint3D CombineTwoNormalsTo(const CPoint3D& pt1, double coef1, const CPoint3D& pt2, double coef2);
}