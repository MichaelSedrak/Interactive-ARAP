#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <iostream>

class Mesh {
public:
    Mesh() {}
    
	void clear(); 

    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getVertices() const;
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getFaces() const;

    void setVertices(const Eigen::MatrixXd& v);

    const int getVertexCount() const;
    const int getFacesCount() const;

    void setVertexAtIndex(int i, double x, double y, double z);

	bool loadMesh(const std::string& filename);
    bool writeMesh(const std::string& filename); 
    
    void verboseOutput();

private:
    // vertices of the mesh
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_vertices;

    // faces of the mesh
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m_faces;
};

