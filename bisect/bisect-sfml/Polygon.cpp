#include "Polygon.h"
#include "geometry.h"
#include <assert.h>

#include <iostream>


Polygon::Polygon() {}


std::vector<Polygon> Polygon::Bisect(float2 p0, float2 p, float2 dir)
{
	int i = GetIndex(p0);
	assert(i != -1);
	return this->BisectInternal(i, p, dir);
}


std::vector<Polygon> Polygon::BisectInternal(int index, float2 p, float2 dir)
{
	int n = (int)this->vertices.size();

	float2 intersection;

	bool found_intersection = false;

	for (int i = 1; i < n; i++)
	{
		if (ray_intersection(p, dir,
			this->vertices[(index + i) % n],
			this->vertices[(index + i + 1) % n],
			intersection))
		{
			found_intersection = true;
			Polygon p1;
			p1.vertices.push_back(this->vertices[index]);
			p1.vertices.push_back(p);
			p1.vertices.push_back(intersection);
			for (int j = (index + i + 1) % n; j != index; j = (j + 1) % n)
			{
				p1.vertices.push_back(this->vertices[j]);
			}

			Polygon p2;
			p2.vertices.push_back(this->vertices[(index + i) % n]);
			p2.vertices.push_back(intersection);
			p2.vertices.push_back(p);
			for (int j = (index + 1) % n; j != (index + i) % n; j = (j + 1) % n)
			{
				p2.vertices.push_back(this->vertices[j]);
			}

			return std::vector<Polygon> { p1, p2 };
		}
	}

	return std::vector<Polygon>();
}


void Polygon::Shrink(float dis)
{

}


void Polygon::Expand(float dis)
{

}

float Polygon::Radius() const
{
	int n = (int)this->vertices.size();
	float x = 0.f, y = 0.f;
	for (auto& v : this->vertices)
	{
		x += v.x;
		y += v.y;
	}
	
	float2 center(x / n, y / n);
	float minDistance = -1.0f;
	for (auto& v : this->vertices)
	{
		float dist = linalg::distance(v, center);
		if (minDistance > dist || minDistance < 0)
			minDistance = dist;
	}

	return minDistance;
}


int Polygon::GetIndex(float2 p)
{
	for (int i = 0; i < this->vertices.size(); ++i)
	{
		if (this->vertices[i] == p) // close to?
		{
			return i;
		}
	}

	return -1;
}