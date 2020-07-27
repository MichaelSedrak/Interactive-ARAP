#pragma once

#include <mesh.h>
#include <solver.h>
#include <demosolver.h>

namespace InteractiveARAP
{
	struct VertexIndexData
	{
		int index;
		double x;
		double y;
		double z;

		bool operator<(const VertexIndexData& a) const
		{
			return index < a.index;
		}
	};

	class NativeInterface
	{
	private:
		std::vector<Mesh> meshes;
		Solver *arapEngine;

		//arap::demo::DemoArapSolver demoArapEngine;
		/*Eigen::VectorXi indices;
		Eigen::MatrixXd positions;
		int meshIDX;*/
		
	public:
		~NativeInterface();
		Mesh deformedMesh;
		NativeInterface();
		void LoadAllMeshes();
		void SetBaseMesh(int index);
		void SetConstraint(int size, int* constraints);
		void SetPositions(int size, float* constraints);
		void Process();
		const double* getMeshVertices(int index) const;
		const int* getMeshIndices(int index) const;
		int getMeshVerticesSize(int index);
		int getMeshIndicesSize(int index);
	};

	class ARAPWrapper
	{
	private:
		static NativeInterface* engine;

	public:
		ARAPWrapper();
		~ARAPWrapper();
		static NativeInterface* GetInstance();
		static void Destroy();
	};


#define EXPORT_API __declspec(dllexport)
	extern "C"
	{
		EXPORT_API void Initialize()
		{
			ARAPWrapper::GetInstance();
			ARAPWrapper::GetInstance()->LoadAllMeshes();
		}

		EXPORT_API const double* GetVertices(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshVertices(modelIndex);
		}

		EXPORT_API const int* GetIndices(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshIndices(modelIndex);
		}

		EXPORT_API int GetIndicesSize(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshIndicesSize(modelIndex);
		}

		EXPORT_API int GetVerticesSize(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshVerticesSize(modelIndex);
		}

		EXPORT_API const double* GetDeformedVertices()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getVertices().data();
		}

		EXPORT_API const int* GetDeformedIndices()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getFaces().data();
		}

		EXPORT_API int GetDeformedIndicesSize()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getFacesCount();
		}

		EXPORT_API int GetDeformedVerticesSize()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getVertexCount();
		}

		EXPORT_API void SetConstraint(int size, int* constraints)
		{
			ARAPWrapper::GetInstance()->SetConstraint(size,constraints);
		}
		
		EXPORT_API void SetPositions(int size, float* constraints)
		{
			ARAPWrapper::GetInstance()->SetPositions(size,constraints);
		}

		EXPORT_API void Process()
		{
			ARAPWrapper::GetInstance()->Process();
		}

		EXPORT_API void SetBaseMesh(int meshIdx)
		{
			ARAPWrapper::GetInstance()->SetBaseMesh(meshIdx);
		}
	}

}
