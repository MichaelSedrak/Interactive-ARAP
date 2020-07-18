#include <iostream>
#include <fstream>

#include "../include/mesh.h"
#include "../include/solver.h"

// g++ -I../libs/eigen-3.3.7/ mesh.cpp solver.cpp test_main_solver.cpp -o test 

int main()
{
	//load mesh
	const std::string filenameSource = std::string("../meshes/armadillo_1k.off");
	Mesh testMesh;

	if (!testMesh.loadMesh(filenameSource)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameSource << std::endl;
		return -1;
	}

	//init ARAP module
	Solver arapSolver = Solver(testMesh.getVertices(), testMesh.getFaces());

    // TODO set constraints
	// arapSolver.SetConstraint()

	// Solve ARAP
	arapSolver.Solve();

	testMesh.setVertices(arapSolver.GetTransformedVertices());
    // Only uncomment if you enjoy scrolling
	// testMesh.verboseOutput();
	testMesh.writeMesh("output.off");

	return 0;
}
