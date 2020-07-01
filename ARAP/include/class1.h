#pragma once

#include <vector>

class dummyClass
{
public:
	dummyClass();
	~dummyClass();
	std::vector<float>& GetVertices();
	std::vector<int>& GetIndices();

private:
	std::vector<float> vertices;
	std::vector<int> indices;

};

