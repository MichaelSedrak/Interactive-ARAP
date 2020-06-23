#include "../include/mesh.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

const unsigned int Vertex::getVertexIdx() const{
    return idx;
}

const Eigen::Vector4f Vertex::getvertexPosition() const{
    return position;
}

void Triangle::calculateFaceNormal(Eigen::Vector4f vert0, Eigen::Vector4f vert1, Eigen::Vector4f vert2){
    Eigen::Vector4f A;
    Eigen::Vector4f B;
    A = vert1 - vert0;
    B = vert2 - vert0; 

    // Nx = Ay * Bz - Az * By
    // Ny = Az * Bx - Ax * Bz
    // Nz = Ax * By - Ay * Bx

    normal.x() = (A.y() * B.z()) - (A.z() * B.y());
    normal.y() = (A.z() * B.x()) - (A.x() * B.z());
    normal.z() = (A.x() * B.y()) - (A.y() * B.x());
    normal.w() = 1.0f;
}

const std::vector<unsigned int> Triangle::getTriangleIdx() const{
    std::vector<unsigned int> triIdx;
    triIdx.push_back(idx0);
    triIdx.push_back(idx1);
    triIdx.push_back(idx2);

    return triIdx;
}

const Eigen::Vector4f Triangle::getNormal() const{
    return normal;
}

void Mesh::clear(){
	m_vertices.clear();
	m_triangles.clear();
} 

// void Mesh::addVertex(Vertex& vertex){
//     vertex.idx = (unsigned int)m_vertices.size();
// 	m_vertices.push_back(vertex);
// } 
// 
// unsigned int Mesh::addFace(unsigned int idx0, unsigned int idx1, unsigned int idx2){
//     Triangle t;
//     t.triangleIdx.x() = idx0;
//     t.triangleIdx.y() = idx1;
//     t.triangleIdx.z() = idx2;
// 	m_triangles.push_back(t);
// } 

const std::vector<Vertex>& Mesh::getVertices() const{
    return m_vertices;
} 

const std::vector<Triangle>& Mesh::getTriangles() const{
    return m_triangles;
} 

void Mesh::transform(const Eigen::Matrix4f& transformation){
	for (Vertex& v : m_vertices) {
		v.position = transformation * v.position;
	}
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
			Vertex v;
			file >> v.position.x() >> v.position.y() >> v.position.z();
			v.position.w() = 1.0f;
            v.idx = i;
			m_vertices.push_back(v);
		}

        // Read faces (i.e. triangles).
    	for (unsigned int i = 0; i < numP; i++) {
    		unsigned int num_vs;
    		file >> num_vs;
            if(num_vs != 3) {
                std::cout << "Only triangles allowed!" << std::endl;
                return false;
            }
    		Triangle t;
    		file >> t.idx0 >> t.idx1 >> t.idx2;
            t.calculateFaceNormal(m_vertices[t.idx0].position, m_vertices[t.idx1].position, m_vertices[t.idx2].position);
    		m_triangles.push_back(t);
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
	outFile << m_vertices.size() << " " << m_triangles.size() << " 0" << std::endl;

	// Save vertices.
	for (unsigned int i = 0; i < m_vertices.size(); i++) {
		const auto& vertex = m_vertices[i];
		if (vertex.position.allFinite())
			outFile << vertex.position.x() << " " << vertex.position.y() << " " << vertex.position.z() << std::endl;
		else
			outFile << "0.0 0.0 0.0" << std::endl;
	}

	// Save faces.
	for (unsigned int i = 0; i < m_triangles.size(); i++) {
		outFile << "3 " << m_triangles[i].idx0 << " " << m_triangles[i].idx1 << " " << m_triangles[i].idx2 << std::endl;
	}

	// Close file.
	outFile.close();
	return true;
} 

void Mesh::verboseOutput(){
    std::cout << "OFF" << std::endl;
    std::cout << m_vertices.size() << " " << m_triangles.size() << " " << "0" << std::endl;
	for (unsigned int i = 0; i < m_vertices.size(); i++) {
		const auto& vertex = m_vertices[i];
		if (vertex.position.allFinite())
			std::cout << vertex.position.x() << " " << vertex.position.y() << " " << vertex.position.z() << std::endl;
		else
			std::cout << "0.0 0.0 0.0" << std::endl;
	}

	for (unsigned int i = 0; i < m_triangles.size(); i++) {
        std::cout << "3 " << m_triangles[i].idx0 << " " << m_triangles[i].idx1 << " " << m_triangles[i].idx2 << std::endl;
	}
}
