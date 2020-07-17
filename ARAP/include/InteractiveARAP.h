#pragma once

#include <class1.h>
#include <mesh.h>

namespace InteractiveARAP
{
	class NativeInterface
	{
	private:
		Mesh meshIO;
		std::vector<std::vector<float>> meshesVertices;
		std::vector<std::vector<int>> meshesIndices;

	public:
		NativeInterface();
		void LoadAllMeshes();
		std::vector<float>& getMeshVertices(int index);
		std::vector<int>& getMeshIndices(int index);
	};

	class ARAP
	{
	private:
		//static dummyClass* engine;
		static NativeInterface* engine;

	public:
		ARAP();
		~ARAP();
		void Process();
		//static dummyClass* GetInstance();
		static NativeInterface* GetInstance();
		static void Destroy();

	};


#define EXPORT_API __declspec(dllexport)
	extern "C"
	{
		EXPORT_API void Initialize()
		{
			ARAP::GetInstance();
			ARAP::GetInstance()->LoadAllMeshes();
		}

		EXPORT_API float* GetVertices(int modelIndex)
		{
			return ARAP::GetInstance()->getMeshVertices(modelIndex).data();
		}

		EXPORT_API int* GetIndices(int modelIndex)
		{
			return ARAP::GetInstance()->getMeshIndices(modelIndex).data();
		}

		EXPORT_API int GetIndicesSize(int modelIndex)
		{
			return ARAP::GetInstance()->getMeshIndices(modelIndex).size();
		}

		EXPORT_API int GetVerticesSize(int modelIndex)
		{
			return ARAP::GetInstance()->getMeshVertices(modelIndex).size();
		}
	}

}