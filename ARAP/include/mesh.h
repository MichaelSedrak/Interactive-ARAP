#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>

struct Vertex {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Eigen::Vector4f position;
    unsigned int idx;

    const unsigned int getVertexIdx() const;
    const Eigen::Vector4f getvertexPosition() const;
};

struct Triangle {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    unsigned int idx0;
    unsigned int idx1;
    unsigned int idx2;
    Eigen::Vector4f normal;

    void calculateFaceNormal(Eigen::Vector4f vert0, Eigen::Vector4f vert1, Eigen::Vector4f vert2);
    const std::vector<unsigned int> getTriangleIdx() const;
    const Eigen::Vector4f getNormal() const;
};

class Mesh {
public:
    Mesh()
    {
    }
    
	void clear(); 
	void addVertex(Vertex& vertex); 
	unsigned int addFace(unsigned int idx0, unsigned int idx1, unsigned int idx2); 

	const std::vector<Vertex>& getVertices() const; 
	const std::vector<Triangle>& getTriangles() const; 

	void transform(const Eigen::Matrix4f& transformation); 
	bool loadMesh(const std::string& filename);
    bool writeMesh(const std::string& filename); 
    
    void verboseOutput();

private:
	std::vector<Vertex> m_vertices;
	std::vector<Triangle> m_triangles;
};

