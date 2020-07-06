#include <iostream>
#include <fstream>

#include "../include/mesh.h"
#include "../include/arap.h"

// g++
// g++ -I../include mesh.cpp test_main.cpp -o test 

int main()
{
	//load mesh
	const std::string filenameSource = std::string("../meshes/armadillo_1k.off");
    Mesh testMesh;

	if (!testMesh.loadMesh(filenameSource)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameSource << std::endl;
		return -1;
	}

    testMesh.verboseOutput();

	//init ARAP module
	ARAP testArap;

	testArap.SetBaseMesh(testMesh);

	std::vector<std::pair<unsigned int, Eigen::Vector3f>> constraints;
	std::pair<unsigned int, Eigen::Vector3f> c0;
	c0 = std::make_pair(0, Eigen::Vector3f(0, 0, 0));
	constraints.push_back(c0);

	testArap.DeformMesh(constraints, 1);

	Mesh deformedMesh = testArap.GetDeformedMesh();

	deformedMesh.verboseOutput();

    return 0;
}
