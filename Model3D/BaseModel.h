// BaseModel.h: interface for the CBaseModel class.
//
//////////////////////////////////////////////////////////////////////
#pragma once

#include "Point3D.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <Eigen/dense>
using namespace std;
namespace Model3D
{
	class CBaseModel
	{
	public:
		CBaseModel();
		CBaseModel(const string& filename);
		struct CFace;
		CBaseModel(const vector<CPoint3D> & vertList, const vector<CFace>& faceList);
		CBaseModel(const vector<CPoint3D> & vertList);
	public:
		struct CFace
		{
			int verts[3];
			CFace() {}
			CFace(int x, int y, int z)
			{
				verts[0] = x;
				verts[1] = y;
				verts[2] = z;
			}
			int& operator[](int index)
			{
				return verts[index];
			}
			int operator[](int index) const
			{
				return verts[index];
			}
			string Encode2String() const
			{
				char buf[1024];
				sprintf_s(buf, "%0d%0d%0d", verts[0], verts[1], verts[2]);
				return string(buf);
			}
			bool operator<(const CFace& other) const
			{
				return Encode2String() < other.Encode2String();
			}
		};
		vector<CFace> m_Faces;
		vector<CPoint3D> m_Verts;

		vector<CPoint3D> m_NormalsToVerts;
		set<int> m_UselessFaces;
		string m_filename;
		double m_scale;

	protected:
		void ReadMFile(const string& filename);
		void ReadFile(const string& filename);
		void ReadObjFile(const string& filename);
		void ReadOffFile(const string& filename);
	public:
		int AddFaceInteriorPoints(int faceIndex, CPoint3D pt);
		map<CPoint3D, int> AddFaceInteriorPoints(const vector<pair<int, CPoint3D>> &pts);
		void PreprocessVertsAndFacesIntoBaseModel();
		double SignedVolume() const;
		void Flip();
		virtual void LoadModel();
		virtual void RemoveUnreferencedVertices();
		int GetVertexID(const CPoint3D& pt) const;
		void SaveMFile(const string& filename) const;
		void SaveOffFile(const string& filename) const;
		void SaveObjFile(const string& filename) const;
		void SaveScalarFieldObjFile(const vector<double>& vals, const string& filename) const;
		void SaveScalarFieldObjFile(const Eigen::VectorXd& vals, const string& filename) const;
		void SaveScalarFieldObjFile(const vector<double>& vals, double maxV, const string& filename) const;
		void SaveScalarFieldObjFile(const Eigen::VectorXd& vals, double maxV, const string& filename) const;
		static void SaveSegmentsToObj(const vector<pair<CPoint3D, CPoint3D>>& segs, const string& filename);
		static void SavePathToObj(const vector<CPoint3D>& pl, const string& filename);
		void SavePamametrizationObjFile(const vector<pair<double, double> >& uvs, const string& filename) const;
		static vector<double> GetScalarField(string filename);
		static string GetComments(const char* filename);
		virtual void PrintInfo(ostream& out) const;
		inline int GetNumOfVerts() const;
		inline int GetNumOfFaces() const;
		string GetFileShortName() const;
		string GetFileFullName() const;
		inline const CPoint3D& Vert(int vertIndex) const;
		inline const CPoint3D& Normal(int vertIndex) const;
		inline const CFace& Face(int faceIndex) const;
		CPoint3D GetShiftVertex(int indexOfVert) const;
		//CPoint3D ShiftVertex(int indexOfVert, double epsilon) const;
		double GetScale() const { return m_scale; }
		vector<CFace> GetUsefulFaces() const;
		inline Eigen::Vector3d GetCoefficients(int faceID, CPoint3D p) const;
	};

	Eigen::Vector3d CBaseModel::GetCoefficients(int faceID, CPoint3D p) const
	{
		Eigen::Matrix<double, 4, 3> M;
		M << Vert(Face(faceID)[0]).x, Vert(Face(faceID)[1]).x, Vert(Face(faceID)[2]).x,
			Vert(Face(faceID)[0]).y, Vert(Face(faceID)[1]).y, Vert(Face(faceID)[2]).y,
			Vert(Face(faceID)[0]).z, Vert(Face(faceID)[1]).z, Vert(Face(faceID)[2]).z,
			1, 1, 1;
		Eigen::Vector4d b;
		b << p.x, p.y, p.z, 1;
		auto MM = M.transpose() * M;
		auto Mb = M.transpose() * b;
		return MM.inverse() * Mb;
	}
	int CBaseModel::GetNumOfVerts() const
	{
		return (int)m_Verts.size();
	}

	int CBaseModel::GetNumOfFaces() const
	{
		return (int)m_Faces.size();
	}

	const CPoint3D& CBaseModel::Vert(int vertIndex) const
	{
		return m_Verts[vertIndex];
	}

	const CPoint3D& CBaseModel::Normal(int vertIndex) const
	{
		return m_NormalsToVerts[vertIndex];
	}

	const CBaseModel::CFace& CBaseModel::Face(int faceIndex) const
	{
		return m_Faces[faceIndex];
	}
}