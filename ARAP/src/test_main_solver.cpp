#include <iostream>
#include <fstream>

#include "../include/mesh.h"
//#include "../include/solver.h"
#include "../include/InteractiveARAP.h"

// g++ -I../libs/eigen-3.3.7/ mesh.cpp solver.cpp test_main_solver.cpp -o test 

int main()
{
	////load mesh
	//const std::string filenameSource = std::string("../meshes/armadillo_1k.off");
	////const std::string filenameSource = std::string("../meshes/square.off");

	//Mesh testMesh;

	//if (!testMesh.loadMesh(filenameSource)) {
	//	std::cout << "Mesh file wasn't read successfully at location: " << filenameSource << std::endl;
	//	return -1;
	//}

	////init ARAP module
	//Solver arapSolver = Solver(testMesh.getVertices(), testMesh.getFaces(),1);

	//for (int i = 0; i < 15; i++)
	//{
	//	arapSolver.Solve();
	//	testMesh.setVertices(arapSolver.GetTransformedVertices());
	//	testMesh.writeMesh("output" + std::to_string(i) + ".off");
	//}
	//return 0;
 //   // TODO set constraints
	//// arapSolver.SetConstraint()

	//// Solve ARAP
	//arapSolver.Solve();

	//testMesh.setVertices(arapSolver.GetTransformedVertices());
 //   // Only uncomment if you enjoy scrolling
	////testMesh.verboseOutput();
	//testMesh.writeMesh("output.off");

	//return 0;

	InteractiveARAP::NativeInterface engine;
	engine.LoadAllMeshes();
	engine.SetBaseMesh(0);
	float positions[32] = {	109.0f,-0.339686f,-0.541938f,0.0259783f,
							328.0f,-0.327594f,-0.56698f,-0.0254948f,
							518.0f,-0.299954f,-0.542617f,0.0895114f,
							316.0f,-0.309019f,-0.568842f,0.0665823f,
							740.0f,-0.239034f,-0.493955f,0.121275f,
							216.0f,-0.221835f,-0.543452f,0.146886f,
							206.0f,-0.242167f,-0.527808f,-0.010471f,
							855.0f,-0.3443436f,-0.3025245f,-0.4035977f};
	engine.SetPositions(32,positions);
	engine.Process();
	return 0;

}

