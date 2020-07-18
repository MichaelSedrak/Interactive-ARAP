#pragma once

#include <mesh.h>
#include <arap.h>

namespace InteractiveARAP
{
	class NativeInterface
	{
	private:
		std::vector<Mesh> meshes;
		ARAP arapEngine;
		
	public:
		Mesh deformedMesh;
		NativeInterface();
		void LoadAllMeshes();
		void Process(int size, float* constraints, unsigned int nIter = 10);
		void SetBaseMesh(int index);
		const std::vector<float>& getMeshVertices(int index);
		const std::vector<int>& getMeshIndices(int index);
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

		EXPORT_API const float* GetVertices(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshVertices(modelIndex).data();
		}

		EXPORT_API const int* GetIndices(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshIndices(modelIndex).data();
		}

		EXPORT_API int GetIndicesSize(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshIndices(modelIndex).size();
		}

		EXPORT_API int GetVerticesSize(int modelIndex)
		{
			return ARAPWrapper::GetInstance()->getMeshVertices(modelIndex).size();
		}

		EXPORT_API const float* GetDeformedVertices()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getVertices().data();
		}

		EXPORT_API const int* GetDeformedIndices()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getTriangles().data();
		}

		EXPORT_API int GetDeformedIndicesSize()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getTriangles().size();
		}

		EXPORT_API int GetDeformedVerticesSize()
		{
			return ARAPWrapper::GetInstance()->deformedMesh.getVertices().size();
		}

		EXPORT_API void Process(int size, float* constraints, unsigned int nIter = 10)
		{
			ARAPWrapper::GetInstance()->Process(size,constraints,nIter);
		}

		EXPORT_API void SetBaseMesh(int meshIdx)
		{
			ARAPWrapper::GetInstance()->SetBaseMesh(meshIdx);
		}
	}

}