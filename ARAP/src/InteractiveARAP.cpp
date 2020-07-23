#include <InteractiveARAP.h>

namespace InteractiveARAP
{
	NativeInterface* ARAPWrapper::engine;

	ARAPWrapper::ARAPWrapper()
	{
	}

	ARAPWrapper::~ARAPWrapper()
	{	
	}

	NativeInterface* ARAPWrapper::GetInstance()
	{
		if (!engine)
			engine = new NativeInterface();
		return engine;
	}

	void ARAPWrapper::Destroy()
	{
		if (engine)
			delete engine;
		engine = nullptr;
	}

	NativeInterface::NativeInterface()
	{
		meshes.resize(10);
	}
	
	void NativeInterface::LoadAllMeshes()
	{
		meshes[0].loadMesh("Assets/meshes/armadillo_1k.off");
		meshes[1].loadMesh("Assets/meshes/bar1.off");
		meshes[2].loadMesh("Assets/meshes/bar2.off");
		meshes[3].loadMesh("Assets/meshes/bar3.off");
		meshes[4].loadMesh("Assets/meshes/cactus_highres.off");
		meshes[5].loadMesh("Assets/meshes/cactus_small.off");
		meshes[6].loadMesh("Assets/meshes/cylinder_small.off");
		meshes[7].loadMesh("Assets/meshes/dino.off");
		meshes[8].loadMesh("Assets/meshes/square_21.off");
		meshes[9].loadMesh("Assets/meshes/square_21_spikes.off");
	}

	void NativeInterface::SetConstraint(int size, int* rawconstraints)
	{
		for (int i = 0; i < size; i++)
			arapEngine.SetConstraint(rawconstraints[i],true);
	}

	void NativeInterface::SetPositions(int size, float* rawconstraints)
	{
		for (int i = 0; i < size; i += 4)
		{
			unsigned int vertexIndex = (unsigned int)rawconstraints[i];
			Eigen::Vector3d vertexConstraint = Eigen::Vector3d((double)rawconstraints[i + 1], (double)rawconstraints[i + 2], (double)rawconstraints[i + 3]);
			arapEngine.SetPosition(vertexIndex,vertexConstraint);
		}
	}

	void NativeInterface::Process()
	{
		arapEngine.Solve();
		deformedMesh.setVertices(arapEngine.GetTransformedVertices());
	}

	void NativeInterface::SetBaseMesh(int index)
	{
		arapEngine = Solver(meshes[index].getVertices(), meshes[index].getFaces());
		deformedMesh = meshes[index];
	}
	
	const double* NativeInterface::getMeshVertices(int index) const
	{
		return meshes[index].getVertices().data();
	}

	const int* NativeInterface::getMeshIndices(int index) const
	{
		return meshes[index].getFaces().data();
	}

	int NativeInterface::getMeshVerticesSize(int index)
	{
		return meshes[index].getVertexCount();
	}

	int NativeInterface::getMeshIndicesSize(int index)
	{
		return meshes[index].getFacesCount();
	}

}