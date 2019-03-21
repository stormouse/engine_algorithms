#include "Scene.h"
#include "geometry.h"


Scene::Scene()
{
	Initialize();
}


void Scene::Render(sf::RenderWindow& window)
{
	if (lines.size() > 0)
	{
		window.draw(&lines[0], lines.size(), sf::Lines);
	}
}

void Scene::Initialize()
{
	Polygon firstPolygon;
	firstPolygon.vertices.push_back(float2(280, 190));
	firstPolygon.vertices.push_back(float2(600, 110));
	firstPolygon.vertices.push_back(float2(760, 300));
	firstPolygon.vertices.push_back(float2(580, 580));
	firstPolygon.vertices.push_back(float2(250, 440));

	this->polygons.push_back(firstPolygon);

	for (int i = 0; i < 5; i++)
	{
		int m = polygons.size();
		for (int j = 0; j < m; j++)
		{
			BisectByLongest(j);
		}
	}

	for (auto polygon : this->polygons)
	{
		int n = (int)polygon.vertices.size();
		for (int i = 0; i < n; i++)
		{
			int j = (i + 1) % n;
			lines.push_back(sf::Vertex(sf::Vector2f(polygon.vertices[i].x, polygon.vertices[i].y), sf::Color::Black));
			lines.push_back(sf::Vertex(sf::Vector2f(polygon.vertices[j].x, polygon.vertices[j].y), sf::Color::Black));
		}
	}
}

void Scene::BisectByLongest(int index)
{
	auto polygon = this->polygons[index];
	
	if (polygon.Radius() < 50.0f) 
		return;
	
	int n = (int)polygon.vertices.size();
	int maxIndex = n - 1;
	float maxLength2 = linalg::distance2(polygon.vertices[n - 1], polygon.vertices[0]);
	for (int i = 0; i < n - 1; i++)
	{
		float l2 = linalg::distance2(polygon.vertices[i], polygon.vertices[i + 1]);
		if (l2 > maxLength2)
		{
			maxLength2 = l2;
			maxIndex = i;
		}
	}

	float2 direction = inner_normal(
		polygon.vertices[maxIndex],
		polygon.vertices[(maxIndex + 1) % n],
		polygon.vertices[(maxIndex + 2) % n]);
	
	float2 p0 = polygon.vertices[maxIndex];
	float2 p1 = polygon.vertices[(maxIndex + 1) % n];
	std::vector<Polygon> result = polygon.Bisect(p0, (p0 + p1) * 0.5f, direction);

	this->polygons.erase(this->polygons.begin() + index);

	for (auto poly : result)
		this->polygons.push_back(poly);
}