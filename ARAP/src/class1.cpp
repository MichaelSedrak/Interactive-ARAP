#include <class1.h>

dummyClass::dummyClass()
{
    vertices.push_back(-1.0f);
    vertices.push_back(-1.0f);
    vertices.push_back(-1.0f);

    vertices.push_back(1.0f);
    vertices.push_back(-1.0f);
    vertices.push_back(-1.0f);

    vertices.push_back(1.0f);
    vertices.push_back(1.0f);
    vertices.push_back(-1.0f);
        
    vertices.push_back(-1.0f);
    vertices.push_back(1.0f);
    vertices.push_back(-1.0f);
    
    vertices.push_back(-1.0f);
    vertices.push_back(-1.0f);
    vertices.push_back(1.0f);

    vertices.push_back(1.0f);
    vertices.push_back(-1.0f);
    vertices.push_back(1.0f);
    
    vertices.push_back(1.0f);
    vertices.push_back(1.0f);
    vertices.push_back(1.0f);
    
    vertices.push_back(-1.0f);
    vertices.push_back(1.0f);
    vertices.push_back(1.0f);

    indices = { 0, 1, 3, 3, 1, 2,
                1, 5, 2, 2, 5, 6,
                5, 4, 6, 6, 4, 7,
                4, 0, 7, 7, 0, 3,
                3, 2, 7, 7, 2, 6,
                4, 5, 0, 0, 5, 1 };
}

dummyClass::~dummyClass()
{
}

std::vector<float>& dummyClass::GetVertices()
{
	return vertices;
}

std::vector<int>& dummyClass::GetIndices()
{
	return indices;
}


