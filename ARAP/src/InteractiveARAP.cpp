#include <InteractiveARAP.h>

namespace InteractiveARAP
{
	//dummyClass* ARAP::engine;
	NativeInterface* ARAP::engine;

	ARAP::ARAP()
	{
	}

	ARAP::~ARAP()
	{	
	}

	void ARAP::Process()
	{
		
	}

	NativeInterface* ARAP::GetInstance()
	{
		if (!engine)
			engine = new NativeInterface();
		return engine;
	}

	void ARAP::Destroy()
	{
		if (engine)
			delete engine;
		engine = nullptr;
	}

	NativeInterface::NativeInterface()
	{
	}
	
	void NativeInterface::LoadAllMeshes()
	{
		meshIO.loadMesh("Assets/meshes/armadillo_1k.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/bar1.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/bar2.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/bar3.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/cactus_highres.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/cactus_small.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/cylinder_small.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/dino.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/square_21.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

		meshIO.loadMesh("Assets/meshes/square_21_spikes.off");
		meshesVertices.push_back(meshIO.getVertices());
		meshesIndices.push_back(meshIO.getTriangles());

	}
	
	std::vector<float>& NativeInterface::getMeshVertices(int index)
	{
		return meshesVertices[index];
	}

	std::vector<int>& NativeInterface::getMeshIndices(int index)
	{
		return meshesIndices[index];
	}
}