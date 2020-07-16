#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>

class Mesh {
public:
    Mesh() {}
    
	void clear(); 

    const Eigen::MatrixXd getVertices() const;
    const Eigen::MatrixXd getFaces() const;

    const int getVertexCount() const;
    const int getFacesCount() const;

    void setVertexAtIndex(int i, double x, double y, double z);

	bool loadMesh(const std::string& filename);
    bool writeMesh(const std::string& filename); 
    
    void verboseOutput();

private:
    // vertices of the mesh
    Eigen::MatrixXd m_vertices;

    // faces of the mesh
    Eigen::MatrixXd m_faces;
};

