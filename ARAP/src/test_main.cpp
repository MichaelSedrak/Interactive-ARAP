#include <iostream>
#include <fstream>

#include "../include/mesh.h"

// g++
// g++ -I../include mesh.cpp test_main.cpp -o test 

int main()
{
	const std::string filenameSource = std::string("../meshes/armadillo_1k.off");
    Mesh testMesh;

	if (!testMesh.loadMesh(filenameSource)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameSource << std::endl;
		return -1;
	}

    testMesh.verboseOutput();

    return 0;
}
