#pragma once

#include "linalg.h"
#include <vector>

using float2 = linalg::float2;

class Polygon
{
public:
	Polygon();

	std::vector<Polygon> Bisect(float2 p0, float2 p, float2 dir);

	void Shrink(float dis);

	void Expand(float dis);

	float Radius() const;

private:
	std::vector<Polygon> BisectInternal(int index, float2 p, float2 dir);

	inline int GetIndex(float2 p);

public:
	std::vector<float2> vertices;
	
};
