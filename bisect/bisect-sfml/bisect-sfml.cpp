#include <SFML/Graphics.hpp>
#include <iostream>
#include "Scene.h"
#include "geometry.h"

int main()
{
	
	sf::RenderWindow window(sf::VideoMode(1024, 768), "Bisect Polygons");

	Scene scene;

	float2 intersection;
	bool intersect = ray_intersection( float2( 413.f, 156.f ), float2( 0.242f, 0.970f ), float2( 520.f, 245.f ), float2( 280.f, 190.f ), intersection);



	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		scene.Step();

		window.clear(sf::Color::White);
		
		scene.Render(window);

		window.display();

	}
}