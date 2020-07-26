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
		meshes[0].loadMesh("../meshes/armadillo_1k.off");
		meshes[1].loadMesh("../meshes/bar1.off");
		meshes[2].loadMesh("../meshes/bar2.off");
		meshes[3].loadMesh("../meshes/bar3.off");
		meshes[4].loadMesh("../meshes/cactus_highres.off");
		meshes[5].loadMesh("../meshes/cactus_small.off");
		meshes[6].loadMesh("../meshes/cylinder_small.off");
		meshes[7].loadMesh("../meshes/dino.off");
		meshes[8].loadMesh("../meshes/square_21.off");
		meshes[9].loadMesh("../meshes/square_21_spikes.off");
	}

	void NativeInterface::SetConstraint(int size, int* rawconstraints)
	{
		/*for (int i = 0; i < size; i++)
			arapEngine.SetConstraint(rawconstraints[i],true);*/
	}

	void NativeInterface::SetPositions(int size, float* rawconstraints)
	{
		/*for (int i = 0; i < size; i += 4)
		{
			unsigned int vertexIndex = (unsigned int)rawconstraints[i];
			arapEngine.SetConstraint(vertexIndex, true);
			Eigen::Vector3d vertexConstraint = Eigen::Vector3d((double)rawconstraints[i + 1], (double)rawconstraints[i + 2], (double)rawconstraints[i + 3]);
			arapEngine.SetPosition(vertexIndex,vertexConstraint);
		}*/
		/*indices.resize(0, 0);
		positions.resize(0, 0);*/

		indices.resize(size / 4, 1);
		positions.resize(size / 4, 3);
		std::vector<VertexIndexData> constraints;
		for (int i = 0; i < size; i += 4)
		{
			VertexIndexData v;
			int vertexIndex = (int)rawconstraints[i];
			v.index = vertexIndex;
			v.x = (double)rawconstraints[i + 1];
			v.y = (double)rawconstraints[i + 2];
			v.z = (double)rawconstraints[i + 3];
			constraints.push_back(v);
		}
		std::sort(constraints.begin(), constraints.end());
		for (int i = 0; i < constraints.size(); i++)
		{
			indices.row(i) << constraints[i].index;
			positions.row(i) << constraints[i].x, constraints[i].y, constraints[i].z;
		}

	}

	void NativeInterface::Process()
	{
		//arapEngine.Solve();
		//deformedMesh.setVertices(arapEngine.GetTransformedVertices());
		//deformedMesh.writeMesh("Assets/meshes/testoutput.off");
		demoArapEngine.RegisterData(meshes[meshIDX].getVertices(), meshes[meshIDX].getFaces(), indices, 15);
		demoArapEngine.Precompute();
		demoArapEngine.Solve(positions);
		deformedMesh.setVertices(demoArapEngine.GetVertexSolution());
		deformedMesh.writeMesh("../meshes/testoutput.off");
	}

	void NativeInterface::SetBaseMesh(int index)
	{
		//arapEngine = Solver(meshes[index].getVertices(), meshes[index].getFaces(),15);
		meshIDX = index;
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
