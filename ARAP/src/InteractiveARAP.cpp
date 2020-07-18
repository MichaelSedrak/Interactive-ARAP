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
		meshes.resize(11);
	}
	
	void NativeInterface::LoadAllMeshes()
	{
		meshes[0].loadMesh("Assets/meshes/armadillo_1k.off");
		meshes[1].loadMesh("Assets/meshes/bar1.off");
		meshes[2].loadMesh("Assets/meshes/bar2.off");
		meshes[3].loadMesh("Assets/meshes/bar2.off");
		meshes[4].loadMesh("Assets/meshes/bar3.off");
		meshes[5].loadMesh("Assets/meshes/cactus_highres.off");
		meshes[6].loadMesh("Assets/meshes/cactus_small.off");
		meshes[7].loadMesh("Assets/meshes/cylinder_small.off");
		meshes[8].loadMesh("Assets/meshes/dino.off");
		meshes[9].loadMesh("Assets/meshes/square_21.off");
		meshes[10].loadMesh("Assets/meshes/square_21_spikes.off");
	}

	void NativeInterface::Process(int size, float* rawconstraints, unsigned int nIter)
	{
		std::vector<std::pair<unsigned int, Eigen::Vector3f>> constraints;
		for (int i = 0; i < size; i += 4)
		{
			unsigned int vertexIndex = (unsigned int)rawconstraints[i];
			Eigen::Vector3f vertexConstraint = Eigen::Vector3f(rawconstraints[i + 1], rawconstraints[i + 2], rawconstraints[i + 3]);
			constraints.push_back(std::pair<unsigned int, Eigen::Vector3f>(vertexIndex, vertexConstraint));
		}

		arapEngine.DeformMesh(constraints, nIter);
		deformedMesh = arapEngine.GetDeformedMesh();
	}

	void NativeInterface::SetBaseMesh(int index)
	{
		arapEngine.SetBaseMesh(meshes[index]);
	}
	
	const std::vector<float>& NativeInterface::getMeshVertices(int index)
	{
		return meshes[index].getVertices();
	}

	const std::vector<int>& NativeInterface::getMeshIndices(int index)
	{
		return meshes[index].getTriangles();
	}
}