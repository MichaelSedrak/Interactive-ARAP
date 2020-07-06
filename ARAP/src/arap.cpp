#pragma once

#include "arap.h"
#include "Eigen.h"
#include "mesh.h"
#include "ProcrustesAligner.h"

//#define USE_DENSE_SYSTEM_MATRIX

ARAP::ARAP() : m_verticesBaseMesh(nullptr), m_verticesDeformed(nullptr), m_rotations(nullptr)
{
	m_procrustesAligner = ProcrustesAligner();
}

ARAP::~ARAP()
{
	// free memory
	SAFE_DELETE_ARRAY(this->m_verticesBaseMesh);
	SAFE_DELETE_ARRAY(this->m_verticesDeformed);
	SAFE_DELETE_ARRAY(this->m_rotations);
}

void ARAP::SetBaseMesh(const Mesh& baseMesh)
{
	// free memory
	SAFE_DELETE_ARRAY(this->m_verticesBaseMesh);
	SAFE_DELETE_ARRAY(this->m_verticesDeformed);
	SAFE_DELETE_ARRAY(this->m_rotations);

	// set base mesh
	this->m_baseMesh = baseMesh;

	// mesh dimensions
	this->m_nVertices = (unsigned int)this->m_baseMesh.getVertexCount();
	this->m_nTriangles = (unsigned int)this->m_baseMesh.getTriangleCount();

	// allocate memory for vertices
	this->m_verticesBaseMesh = new Eigen::Vector3f[this->m_nVertices];
	this->m_verticesDeformed = new Eigen::Vector3f[this->m_nVertices];
	std::vector<float> vertexList = baseMesh.getVertices();
	for (unsigned int i = 0; i < this->m_nVertices; i++)
	{
		this->m_verticesBaseMesh[i] = Eigen::Vector3f(vertexList[3*i], vertexList[3 * i + 1], vertexList[3 * i + 2]);
	}
	memcpy(this->m_verticesDeformed, this->m_verticesBaseMesh, sizeof(Eigen::Vector3f) * this->m_nVertices);

	// allocate memory for rotations
	this->m_rotations = new Eigen::Matrix3f[this->m_nVertices];

	// init system matrix
	InitSystemMatrix();
}

void ARAP::DeformMesh(const std::vector<std::pair<unsigned int, Eigen::Vector3f>>& constraints, unsigned int nIter = 10)
{
	VERBOSE("Deform mesh ...");
	// reset to unmodified mesh
	memcpy(this->m_verticesDeformed, this->m_verticesBaseMesh, sizeof(Eigen::Vector3f) * this->m_nVertices);

	// apply constraints
	if (constraints.empty()) return;
	for (auto c : constraints)
	{
		this->m_verticesDeformed[c.first] = c.second;
	}

	// run flip-flop optimization
	ProgressBar("Deform mesh", 0.0f);
	for (unsigned int iter = 0; iter < nIter; ++iter)
	{
		// compute rotations
		SolveForRotations();

		// compute optimal vertex positions
		SolveForVertexPositions(constraints);

		ProgressBar("Deform mesh", iter / (nIter - 1.0f));
	}

	VERBOSE("Deform mesh ... DONE!");
}

Mesh ARAP::GetDeformedMesh()
{
	// copy base mesh
	Mesh result = this->m_baseMesh;

	// set deformed vertices
	for (unsigned int i = 0; i < this->m_nVertices; i++)
	{
		result.setVertexAtIndex(i, this->m_verticesDeformed[i].x(), this->m_verticesDeformed[i].y(), this->m_verticesDeformed[i].z());
	}

	return result;
}

void ARAP::InitSystemMatrix()
{
	VERBOSE("Init system matrix ...");

	// allocate memory
#ifdef USE_DENSE_SYSTEM_MATRIX
	this->m_systemMatrix = Eigen::MatrixXf(this->m_nVertices, this->m_nVertices);
	this->m_systemMatrix.setZero();
#else
	this->m_systemMatrixSparse = Eigen::SparseMatrix<float>(this->m_nVertices, this->m_nVertices);
#endif
	this->m_rhsX = Eigen::VectorXf(this->m_nVertices);
	this->m_rhsY = Eigen::VectorXf(this->m_nVertices);
	this->m_rhsZ = Eigen::VectorXf(this->m_nVertices);

	// reset matrix
	this->m_rhsX.setZero();
	this->m_rhsY.setZero();
	this->m_rhsZ.setZero();

	std::vector<int> triangleList = this->m_baseMesh.getTriangles();

	// fill matrix
	for (unsigned int i = 0; i < this->m_nVertices; ++i)
	{
		this->m_systemMatrixSparse.insert(i, i) = 0.0f;
		//for (Mesh::VertexVertexIter vv_it = this->m_baseMesh.vv_iter(VertexHandle(i)); vv_it; ++vv_it)

		for(int k = 0; k < this->m_nTriangles; k++)
		{
			int v1 = triangleList[k * 3];
			int v2 = triangleList[k * 3 + 1];
			int v3 = triangleList[k * 3 + 2];

			if (v1 == i) {
				float w_ij = 1.0f;
#ifdef USE_DENSE_SYSTEM_MATRIX
				//add edge to v2
				this->m_systemMatrix(i, i) += w_ij;
				this->m_systemMatrix(i, v2) = -w_ij;
#else
				this->m_systemMatrixSparse.coeffRef(i, i) += w_ij;
				this->m_systemMatrixSparse.insert(i, v2) = -w_ij;
					
				//add edge to v3
				this->m_systemMatrix(i, i) += w_ij;
				this->m_systemMatrix(i, v3) = -w_ij;
#else
				this->m_systemMatrixSparse.coeffRef(i, i) += w_ij;
				this->m_systemMatrixSparse.insert(i, v3) = -w_ij;
#endif
			}
			else if (v2 == i) {
				float w_ij = 1.0f;
#ifdef USE_DENSE_SYSTEM_MATRIX
				//add edge to v1
				this->m_systemMatrix(i, i) += w_ij;
				this->m_systemMatrix(i, v1) = -w_ij;
#else
				this->m_systemMatrixSparse.coeffRef(i, i) += w_ij;
				this->m_systemMatrixSparse.insert(i, v1) = -w_ij;

				//add edge to v3
				this->m_systemMatrix(i, i) += w_ij;
				this->m_systemMatrix(i, v3) = -w_ij;
#else
				this->m_systemMatrixSparse.coeffRef(i, i) += w_ij;
				this->m_systemMatrixSparse.insert(i, v3) = -w_ij;
#endif
			}
			else if (v3 == i) {
				float w_ij = 1.0f;
#ifdef USE_DENSE_SYSTEM_MATRIX
				//add edge to v1
				this->m_systemMatrix(i, i) += w_ij;
				this->m_systemMatrix(i, v1) = -w_ij;
#else
				this->m_systemMatrixSparse.coeffRef(i, i) += w_ij;
				this->m_systemMatrixSparse.insert(i, v1) = -w_ij;

				//add edge to v3
				this->m_systemMatrix(i, i) += w_ij;
				this->m_systemMatrix(i, v2) = -w_ij;
#else
				this->m_systemMatrixSparse.coeffRef(i, i) += w_ij;
				this->m_systemMatrixSparse.insert(i, v2) = -w_ij;
#endif
			}
		}
	}

	VERBOSE("Init system matrix ... DONE!");
}

void ARAP::ComputeRightHandSide()
{
	VERBOSE("Compute right hand side ...");

	// reset rhs
	this->m_rhsX.setZero();
	this->m_rhsY.setZero();
	this->m_rhsZ.setZero();

#pragma omp parallel for
	for (int i = 0; i < (int)this->m_nVertices; ++i)
	{
		Eigen::Vector3f sum(0.0f, 0.0f, 0.0f);

		std::vector<int> triangleList = this->m_baseMesh.getTriangles();

		//for (Mesh::VertexVertexIter vv_it = this->m_baseMesh.vv_iter(VertexHandle(i)); vv_it; ++vv_it)
		for (int k = 0; k < this->m_nTriangles; k++)
		{
			int v1 = triangleList[k * 3];
			int v2 = triangleList[k * 3 + 1];
			int v3 = triangleList[k * 3 + 2];

			float w_ij = 1.0f;
			if (v1 == i) {
				sum += w_ij / 2.0f * (this->m_rotations[i] + this->m_rotations[v2]) * (this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v2]);
				sum += w_ij / 2.0f * (this->m_rotations[i] + this->m_rotations[v3]) * (this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v3]);
			}
			else if (v2 == i) {
				sum += w_ij / 2.0f * (this->m_rotations[i] + this->m_rotations[v1]) * (this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v1]);
				sum += w_ij / 2.0f * (this->m_rotations[i] + this->m_rotations[v3]) * (this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v3]);
			}
			else if (v3 == i) {
				sum += w_ij / 2.0f * (this->m_rotations[i] + this->m_rotations[v1]) * (this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v1]);
				sum += w_ij / 2.0f * (this->m_rotations[i] + this->m_rotations[v2]) * (this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v2]);
			}
		}
		this->m_rhsX(i) = sum.x();
		this->m_rhsY(i) = sum.y();
		this->m_rhsZ(i) = sum.z();
	}

	VERBOSE("Compute right hand side ... DONE!");
}

void ARAP::SolveForRotations()
{
	VERBOSE("Solve for rotations ...");

#pragma omp parallel for
	for (int i = 0; i < (int)this->m_nVertices; ++i)
	{
		std::vector<Eigen::Vector3f> basePositions;
		std::vector<Eigen::Vector3f> currentPositions;

		std::vector<int> triangleList = this->m_baseMesh.getTriangles();

		//for (Mesh::VertexVertexIter vv_it = this->m_baseMesh.vv_iter(VertexHandle(i)); vv_it; ++vv_it)
		for (int k = 0; k < this->m_nTriangles; k++)
		{
			int v1 = triangleList[k * 3];
			int v2 = triangleList[k * 3 + 1];
			int v3 = triangleList[k * 3 + 2];

			if (v1 == i) {
				basePositions.push_back(this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v2]);
				currentPositions.push_back(this->m_verticesDeformed[i] - this->m_verticesDeformed[v2]);

				basePositions.push_back(this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v3]);
				currentPositions.push_back(this->m_verticesDeformed[i] - this->m_verticesDeformed[v3]);
			}
			else if (v2 == i) {
				basePositions.push_back(this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v1]);
				currentPositions.push_back(this->m_verticesDeformed[i] - this->m_verticesDeformed[v1]);

				basePositions.push_back(this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v3]);
				currentPositions.push_back(this->m_verticesDeformed[i] - this->m_verticesDeformed[v3]);
			}
			else if (v3 == i) {
				basePositions.push_back(this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v1]);
				currentPositions.push_back(this->m_verticesDeformed[i] - this->m_verticesDeformed[v1]);

				basePositions.push_back(this->m_verticesBaseMesh[i] - this->m_verticesBaseMesh[v2]);
				currentPositions.push_back(this->m_verticesDeformed[i] - this->m_verticesDeformed[v2]);
			}
		}

		this->m_rotations[i] = m_procrustesAligner.estimatePose(currentPositions, basePositions).block<3, 3>(0, 0);
	}

	VERBOSE("Solve for rotations ... DONE!");
}

void ARAP::SolveForVertexPositions(const std::vector<std::pair<unsigned int, Eigen::Vector3f>>& constraints)
{
	VERBOSE("Solve for rotations ...");

	// compute rhs
	ComputeRightHandSide();

	// solve system
#ifdef USE_DENSE_SYSTEM_MATRIX
	Eigen::MatrixXf systemMatrix = this->m_systemMatrix;
#else
	Eigen::SparseMatrix<float> systemMatrixSparse = this->m_systemMatrixSparse;
#endif

	// adjust system matrix
	for (auto c : constraints)
	{
		// adapt rhs
		for (unsigned int i = 0; i < this->m_nVertices; ++i)
		{
#ifdef USE_DENSE_SYSTEM_MATRIX
			this->m_rhsX(i) -= systemMatrix(i, c.first) * c.second.x();
			this->m_rhsY(i) -= systemMatrix(i, c.first) * c.second.y();
			this->m_rhsZ(i) -= systemMatrix(i, c.first) * c.second.z();
#else
			if (systemMatrixSparse.coeff(i, c.first) != 0.0f)
			{
				this->m_rhsX(i) -= systemMatrixSparse.coeff(i, c.first) * c.second.x();
				this->m_rhsY(i) -= systemMatrixSparse.coeff(i, c.first) * c.second.y();
				this->m_rhsZ(i) -= systemMatrixSparse.coeff(i, c.first) * c.second.z();
			}
#endif
		}
		this->m_rhsX(c.first) = c.second.x();
		this->m_rhsY(c.first) = c.second.y();
		this->m_rhsZ(c.first) = c.second.z();

		// delete row and column
#ifdef USE_DENSE_SYSTEM_MATRIX
		for (unsigned int i = 0; i < this->m_nVertices; ++i) systemMatrix(c.first, i) = systemMatrix(i, c.first) = 0.0f;
		systemMatrix(c.first, c.first) = 1.0f;
#else
		for (unsigned int i = 0; i < this->m_nVertices; ++i) if (systemMatrixSparse.coeff(c.first, i) != 0.0f) systemMatrixSparse.coeffRef(c.first, i) = systemMatrixSparse.coeffRef(i, c.first) = 0.0f;
		systemMatrixSparse.coeffRef(c.first, c.first) = 1.0f;
#endif

	}

#ifdef USE_DENSE_SYSTEM_MATRIX
	static Eigen::JacobiSVD<Eigen::MatrixXf> svd(systemMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
#else
	/*static*/ Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> svd(systemMatrixSparse);
#endif
	Eigen::VectorXf x = svd.solve(this->m_rhsX);
	Eigen::VectorXf y = svd.solve(this->m_rhsY);
	Eigen::VectorXf z = svd.solve(this->m_rhsZ);

	for (unsigned int i = 0; i < this->m_nVertices; ++i)
	{
		this->m_verticesDeformed[i] = Eigen::Vector3f(x(i), y(i), z(i));
	}

	VERBOSE("Solve for rotations ... DONE!");
}

void ARAP::ProgressBar(std::string titel, float progress)
{
	int barWidth = 70;

	std::cout << titel << " [";
	int pos = int(barWidth * progress);
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();

	if (progress == 1.0f)	std::cout << std::endl;
}