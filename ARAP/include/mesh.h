#pragma once

#include <iostream>
#include <vector>

class Mesh {
public:
    Mesh() {}
    
	void clear(); 

    const std::vector<float>& getVertices() const;
    const std::vector<int>& getTriangles() const;

    const int getVertexCount() const; //TODO
    const int getTriangleCount() const; //TODO

    void setVertexAtIndex(int i, float x, float y, float z); //TODO

    void transform(const float& transformation);
	bool loadMesh(const std::string& filename);
    bool writeMesh(const std::string& filename); 
    
    void verboseOutput();

private:
    std::vector<float> m_vertices;
    std::vector<int> m_triangles;
};

