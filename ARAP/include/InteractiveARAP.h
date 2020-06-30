#pragma once

#include <class1.h>

namespace InteractiveARAP
{
	class ARAP
	{
	private:
		static dummyClass* engine;

	public:
		ARAP();
		~ARAP();
		void Process();
		static dummyClass* GetInstance();
		static void Destroy();

	};


#define EXPORT_API __declspec(dllexport)
	extern "C"
	{
		EXPORT_API void Initialize()
		{
			ARAP::GetInstance();
		}

		EXPORT_API float* GetVertices()
		{
			return ARAP::GetInstance()->GetVertices().data();
		}

		EXPORT_API int* GetIndices()
		{
			return ARAP::GetInstance()->GetIndices().data();
		}

		EXPORT_API int GetIndicesSize()
		{
			return ARAP::GetInstance()->GetIndices().size();
		}

		EXPORT_API int GetVerticesSize()
		{
			return ARAP::GetInstance()->GetVertices().size();
		}
	}

}