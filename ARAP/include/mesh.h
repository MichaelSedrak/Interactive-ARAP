#pragma once

#include <iostream>
#include <vector>

class Mesh {
public:
    Mesh() {}
    
	void clear(); 

    std::vector<float>& getVertices();
    std::vector<int>& getTriangles();

    void transform(const float& transformation);
	bool loadMesh(const std::string& filename);
    bool writeMesh(const std::string& filename); 
    
    void verboseOutput();

private:
    std::vector<float> m_vertices;
    std::vector<int> m_triangles;
};

