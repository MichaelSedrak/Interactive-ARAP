#include "../include/mesh.h"
#include <fstream>
#include <iostream>

void Mesh::clear(){
    // Deallocate memory
	m_vertices.resize(0, 0);
	m_faces.resize(0, 0);
} 

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& Mesh::getVertices() const{
    return m_vertices;
} 

const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& Mesh::getFaces() const{
    return m_faces;
} 

void Mesh::setVertices(const Eigen::MatrixXd& v) {
	m_vertices = v;
}

const int Mesh::getVertexCount() const{
	return m_vertices.rows();
}

const int Mesh::getFacesCount() const{
	return m_faces.rows();
}

void Mesh::setVertexAtIndex(int i, double x, double y, double z) {
    m_vertices.row(i) << x, y, z;
}

bool Mesh::loadMesh(const std::string& filename){
    std::cout << ".off file: loading..." << std::endl;
    // Read off file (Important: Only .off files are supported).
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

    m_vertices.resize(numV, 3);
    m_faces.resize(numP, 3);

    if (std::string(string1).compare("OFF") == 0) {
	    // Read vertices.
		for (unsigned int i = 0; i < numV; i++) {
			float x, y, z;
            file >> x >> y >> z;
            m_vertices.row(i) << x, y, z;
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
            m_faces.row(i) << idx0, idx1, idx2;
    	}
	}
	else {
		std::cout << "Incorrect mesh file type." << std::endl;
		return false;
	}
    std::cout << ".off file: loading done." << std::endl;
	return true;
}

bool Mesh::writeMesh(const std::string& filename){
	// Write off file.
	std::ofstream outFile(filename);
	if (!outFile.is_open()) return false;

	// Write header.
	outFile << "OFF" << std::endl;
	outFile << m_vertices.rows() << " " << m_faces.rows() << " 0" << std::endl;

	// Save vertices.
	for (unsigned int i = 0; i < m_vertices.rows(); i++) {
		outFile << m_vertices(i, 0) << " " << m_vertices(i, 1) << " " << m_vertices(i, 2) << std::endl;
	}

	// Save faces.
	for (unsigned int i = 0; i < m_faces.rows(); i++) {
		outFile << "3 " << m_faces(i, 0) << " " << m_faces(i, 1) << " " << m_faces(i, 2) << std::endl;
	}

	// Close file.
	outFile.close();
	return true;
} 

void Mesh::verboseOutput(){
    // Write header.
    std::cout << "OFF" << std::endl;
    std::cout << m_vertices.rows() << " " << m_faces.rows() << " 0" << std::endl;

	// Save vertices.
	for (unsigned int i = 0; i < m_vertices.rows(); i++) {
        std::cout << m_vertices(i, 0) << " " << m_vertices(i, 1) << " " << m_vertices(i, 2) << std::endl;
	}

	// Save faces.
	for (unsigned int i = 0; i < m_faces.rows(); i++) {
        std::cout << "3 " << m_faces(i, 0) << " " << m_faces(i, 1) << " " << m_faces(i, 2) << std::endl;
	}
}
