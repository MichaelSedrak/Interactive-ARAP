#include <InteractiveARAP.h>

namespace InteractiveARAP
{
	dummyClass* ARAP::engine;

	ARAP::ARAP()
	{
	}

	ARAP::~ARAP()
	{	
	}

	void ARAP::Process()
	{
		
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