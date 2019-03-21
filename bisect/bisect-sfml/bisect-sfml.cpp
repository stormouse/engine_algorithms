#include <SFML/Graphics.hpp>
#include <iostream>
#include "Scene.h"

int main()
{
	
	sf::RenderWindow window(sf::VideoMode(1024, 768), "Bisect Polygons");

	Scene scene;

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear(sf::Color::White);
		
		scene.Render(window);

		window.display();
	}
}