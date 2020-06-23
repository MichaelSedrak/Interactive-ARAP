#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include "../include/mesh.h"

// g++ $(pkg-config --cflags eigen3) -I../include mesh.cpp test_main.cpp -o test 

int main()
{
	const std::string filenameSource = std::string("../meshes/bar1.off");
    Mesh testMesh;

	if (!testMesh.loadMesh(filenameSource)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameSource << std::endl;
		return -1;
	}

    //std::cout << testMesh.getVertices() << std::endl;
    //std::cout << testMesh.getTriangles() << std::endl;

    testMesh.verboseOutput();

    return 0;
}
