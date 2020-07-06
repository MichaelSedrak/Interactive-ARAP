#pragma once

#include "Eigen.h"
#include "mesh.h"
#include "ProcrustesAligner.h"

//#define USE_DENSE_SYSTEM_MATRIX

class ARAP
{
public:
	ARAP() : m_verticesBaseMesh(nullptr), m_verticesDeformed(nullptr), m_rotations(nullptr){}

	~ARAP() {}

	void SetBaseMesh(const Mesh& baseMesh);
	void DeformMesh(const std::vector<std::pair<unsigned int, Eigen::Vector3f>>& constraints, unsigned int nIter = 10);
	Mesh GetDeformedMesh();

private:
	void InitSystemMatrix();

	void ComputeRightHandSide();

	void SolveForRotations();

	void SolveForVertexPositions(const std::vector<std::pair<unsigned int, Eigen::Vector3f>>& constraints);

	void ProgressBar(std::string titel, float progress);

	// base mesh
	Mesh m_baseMesh;

	// nVertices
	unsigned int m_nVertices;

	// nTriangles
	unsigned int m_nTriangles;

	// vertices of the base mesh
	Eigen::Vector3f* m_verticesBaseMesh;

	// vertices of the current deformed mesh
	Eigen::Vector3f* m_verticesDeformed;

	// rotations of the vertices
	Eigen::Matrix3f* m_rotations;

	// system matrix
#ifdef USE_DENSE_SYSTEM_MATRIX
	Eigen::MatrixXf m_systemMatrix;
#else
	Eigen::SparseMatrix<float> m_systemMatrixSparse;
#endif

	// right hand side
	Eigen::VectorXf m_rhsX;
	Eigen::VectorXf m_rhsY;
	Eigen::VectorXf m_rhsZ;

	ProcrustesAligner m_procrustesAligner;
};
