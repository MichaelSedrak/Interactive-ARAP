#include <InteractiveARAP.h>
namespace InteractiveARAP
{
	ARAP::ARAP()
	{
		//engine = new dummyClass();
	}

	ARAP::~ARAP()
	{
		//delete engine;
	}

	void ARAP::Process()
	{
		std::vector<float>& nativeVertices = engine->GetVertices();
		std::vector<int>& nativeIndices = engine->GetIndices();

		vertices = gcnew System::Collections::Generic::List<float>(nativeVertices.size());
		indices = gcnew System::Collections::Generic::List<int>(nativeIndices.size());

		for (int i = 0; i < nativeVertices.size(); i++)
			vertices->Add(nativeVertices[i]);

		for (int i = 0; i < nativeIndices.size(); i++)
			indices->Add(nativeIndices[i]);
	}

	dummyClass* ARAP::GetInstance()
	{
		if (!engine)
			engine = new dummyClass();
		return engine;
	}

	void ARAP::Destroy()
	{
		if (engine)
			delete engine;
		engine = nullptr;
	}
}