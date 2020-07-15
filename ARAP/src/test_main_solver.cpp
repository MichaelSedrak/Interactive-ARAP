#include <iostream>
#include <fstream>

#include "../include/mesh.h"
#include "../include/solver.h"

// g++
// g++ -I../include -I../libs/libigl-2.2.0/include -I../libs/eigen-3.3.7/ mesh.cpp arap.cpp test_main_refacto_refactor.cpp -o test

int main()
{
	//load mesh
	const std::string filenameSource = std::string("../meshes/armadillo_1k.off");
	Mesh testMesh;

	if (!testMesh.loadMesh(filenameSource)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameSource << std::endl;
		return -1;
	}

	// testMesh.verboseOutput();

	//init ARAP module
	Solver arapSolver;

	arapSolver.LoadData(testMesh);


	// Precompute 1-Ring-Neigborhood
	arapSolver.PrecomputeNeighbors();

	// Precomputes the weights used in the paper
	arapSolver.PrecomputeCotangentWeights();

	// Compute Cholesky factorization of L
	arapSolver.PrecomputeLaplace-Beltrami();

	/*
	//Setup Constraints
	std::vector<std::pair<unsigned int, Eigen::Vector3f>> constraints;
	std::pair<unsigned int, Eigen::Vector3f> c0;
	//-0.317609 -0.547101 0.0205811 
	c0 = std::make_pair(0, Eigen::Vector3f(0.3176091, -0.547101, 0.0205811));
	constraints.push_back(c0);
	*/

	// Solve ARAP
	arapSolver.Solve();

	//TODO
	//Mesh deformedMesh = arapSolver.SOMEHOW_GET_MESH();
	//deformedMesh.writeMesh("output.off");

	return 0;
}
