#include "../include/mesh.h"
#include <fstream>
#include <iostream>
#include <vector>

void Mesh::clear(){
	m_vertices.clear();
	m_triangles.clear();
} 

const std::vector<float>& Mesh::getVertices() const{
    return m_vertices;
} 

const std::vector<int>& Mesh::getTriangles() const{
    return m_triangles;
} 

void Mesh::transform(const float& transformation){
    // TODO
}

const int Mesh::getVertexCount() const{
	return m_vertices.size()/3;
}
const int Mesh::getTriangleCount() const{
	return m_triangles.size() / 3;
}

void Mesh::setVertexAtIndex(int i, float x, float y, float z) {
	//TODO
}

bool Mesh::loadMesh(const std::string& filename){
    // Read off file (Important: Only .off files are supported).
	m_vertices.clear();
	m_triangles.clear();

	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "Mesh file wasn't read successfully." << std::endl;
		return false;
	}

	// First line should say 'OFF'.
	char string1[5];
	file >> string1;

	// Read header.
	unsigned int numV = 0;
	unsigned int numP = 0;
	unsigned int numE = 0;
	file >> numV >> numP >> numE;

	m_vertices.reserve(numV);
	m_triangles.reserve(numP);

    if (std::string(string1).compare("OFF") == 0) {
	    // Read vertices.
		for (unsigned int i = 0; i < numV; i++) {
			float x, y, z;
            file >> x >> y >> z;
			m_vertices.push_back(x);
			m_vertices.push_back(y);
			m_vertices.push_back(z);
		}

        // Read faces (i.e. triangles).
    	for (unsigned int i = 0; i < numP; i++) {
    		unsigned int num_vs;
    		file >> num_vs;
            if(num_vs != 3) {
                std::cout << "Only triangles allowed!" << std::endl;
                return false;
            }
            int idx0, idx1, idx2;
    		file >> idx0 >> idx1 >> idx2;
    		m_triangles.push_back(idx0);
    		m_triangles.push_back(idx1);
    		m_triangles.push_back(idx2);
    	}
	}
	else {
		std::cout << "Incorrect mesh file type." << std::endl;
		return false;
	}

	return true;
}

bool Mesh::writeMesh(const std::string& filename){
	// Write off file.
	std::ofstream outFile(filename);
	if (!outFile.is_open()) return false;

	// Write header.
	outFile << "OFF" << std::endl;
	outFile << (m_vertices.size() / 3) << " " << (m_triangles.size() / 3) << " 0" << std::endl;

	// Save vertices.
	for (unsigned int i = 0; i < m_vertices.size(); i += 3) {
		outFile << m_vertices[i] << " " << m_vertices[i+1] << " " << m_vertices[i+2] << std::endl;
	}

	// Save faces.
	for (unsigned int i = 0; i < m_triangles.size(); i += 3) {
		outFile << "3 " << m_triangles[i] << " " << m_triangles[i+1] << " " << m_triangles[i+2] << std::endl;
	}

	// Close file.
	outFile.close();
	return true;
} 

void Mesh::verboseOutput(){
    // Write header.
    std::cout << "OFF" << std::endl;
    std::cout << (m_vertices.size() / 3) << " " << (m_triangles.size() / 3) << " 0" << std::endl;

	// Save vertices.
	for (unsigned int i = 0; i < m_vertices.size(); i += 3) {
        std::cout << m_vertices[i] << " " << m_vertices[i+1] << " " << m_vertices[i+2] << std::endl;
	}

	// Save faces.
	for (unsigned int i = 0; i < m_triangles.size(); i += 3) {
        std::cout << "3 " << m_triangles[i] << " " << m_triangles[i+1] << " " << m_triangles[i+2] << std::endl;
	}


}
