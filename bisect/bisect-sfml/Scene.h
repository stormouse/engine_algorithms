#pragma once

#include <vector>
#include <SFML/Graphics.hpp>

#include "Polygon.h"

class Scene
{

public:
	
	Scene();

	void Render(sf::RenderWindow& window);


private:

	void Initialize();

	void BisectByLongest(int index);

private:
	std::vector<Polygon> polygons;
	std::vector<sf::Vertex> lines;

};
