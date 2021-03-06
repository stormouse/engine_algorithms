#include <SFML/Graphics.hpp>

#include "linalg.h"
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <random>
#include <chrono>
#include <memory>
#include <cmath>

#define min(a, b)(((a) < (b)) ? (a) : (b))
#define max(a, b)(((a) > (b)) ? (a) : (b))


class Random 
{
public:
	Random() : Random(0) {}
	Random(unsigned int seed) : generator(seed), distribution(0.0f, 1.0f) {}

	float randomFloat(float min, float max) 
	{
		return distribution(generator) * (max - min) + min;
	}

private:
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution;
};

struct Context
{
	sf::RenderWindow *pWindow;
};

struct ProblemConfig
{
	int numPoints;
	int spaceWidth;
	int spaceHeight;
	int spaceX;
	int spaceY;
};

struct Triangle;
struct Edge;

typedef std::list<Edge>::iterator EdgeNode;
typedef std::list<Triangle>::iterator TriangleNode;
typedef std::list<linalg::float2>::iterator VertexNode;

struct Edge
{
	linalg::float2 a;
	linalg::float2 b;
	TriangleNode pTriangle[2];

	Edge(const linalg::float2& _a, const linalg::float2& _b) : a(_a), b(_b) {}
};

struct Triangle
{
	EdgeNode pEdge[3];
};

class Problem
{
public:
	Problem(Context& context, const ProblemConfig& config) : m_context(context), random(1337) 
	{
		configure(config);
	}

	void start()
	{
		generateRandomVertices();
	}

	void drawGeneratedVertices()
	{
		for (auto& v : vertices)
		{
			sf::CircleShape point(2.5f);
			point.setFillColor(sf::Color::Green);
			point.setPosition(sf::Vector2f(v.x, v.y));
			m_context.pWindow->draw(point);
		}
	}

	std::vector<linalg::float2> getVertices() const
	{
		return vertices;
	}

private:
	void configure(const ProblemConfig& config)
	{
		this->numPoints = config.numPoints;
		this->spaceWidth = config.spaceWidth;
		this->spaceHeight = config.spaceHeight;
		this->spaceX = config.spaceX;
		this->spaceY = config.spaceY;
	}

	void generateRandomVertices()
	{
		for (int i = 0; i < numPoints; i++)
		{
			int q = sqrt(numPoints) + 1;
			float dx = 1.f * spaceWidth / q;
			float dy = 1.f * spaceHeight / q;
			vertices.push_back(linalg::float2(
				spaceX + (i % q) * dx + random.randomFloat(dx * -0.45f, dx * 0.45f),
				spaceY + (i / q) * dy + random.randomFloat(dy * -0.45f, dy * 0.45f)
			));
		}
	}

	void testPrintGeneratedVertices()
	{
		for (auto& v : vertices)
		{
			std::cout << v.x << ", " << v.y << "\n";
		}
	}

	void testDrawGeneratedVertices()
	{
		for (auto& v : vertices)
		{
			sf::CircleShape point(2.5f);
			point.setFillColor(sf::Color::Green);
			point.setPosition(sf::Vector2f(v.x, v.y));
			m_context.pWindow->draw(point);
		}
	}

private:
	Context m_context;

	int numPoints;
	int spaceWidth;
	int spaceHeight;
	int spaceX;
	int spaceY;

	std::vector<linalg::float2> vertices;
	Random random;
};

struct TriangulationResult
{
	int error;
	std::vector<Edge> edges;
};

class DelaunayTriangulation
{

public:
	TriangulationResult Triangulate(const std::vector<linalg::float2>& vertices)
	{
		TriangulationResult result;

		feedVertices(vertices);
		generateBoundTriangle();
		for (const auto& v : vertices)
			addVertex(v);
		
		VertexNode vn = --this->vertices.end();
		linalg::float2 v1 = *vn; vn--;
		linalg::float2 v2 = *vn; vn--;
		linalg::float2 v3 = *vn; vn--;

		for (const auto& e : edges)
		{
			if (e.a == v1 || e.a == v2 || e.a == v3) continue;
			if (e.b == v1 || e.b == v2 || e.b == v3) continue;
			result.edges.push_back(e);
		}

		result.error = 0;

		return result;
	}

	void TriangulateSetup(const std::vector<linalg::float2>& vertices)
	{
		feedVertices(vertices);
		generateBoundTriangle();
		vtxit = this->vertices.begin();
		vtxcount = vertices.size();
		vtxitctr = 0;
	}

	TriangulationResult TriangulateStep()
	{
		TriangulationResult result;
		if (vtxitctr < vtxcount)
		{
			addVertex(*vtxit);
			vtxit++;
			vtxitctr++;
		}

		result.error = 0;
		for (const auto& e : edges)
			result.edges.push_back(e);
		
		return result;
	}

	VertexNode vtxit;
	int vtxcount;
	int vtxitctr;

private:
	void feedVertices(const std::vector<linalg::float2>& _vertices)
	{
		vertices.clear();
		for (auto& v : _vertices)
			vertices.push_back(v);
	}

	void generateBoundTriangle()
	{
		if (vertices.size() < 3) throw "vertices.size() < 3";

		float xmin, xmax, ymin, ymax;
		xmin = xmax = vertices.front().x;
		ymin = ymax = vertices.front().y;
		for (VertexNode it = vertices.begin(); it != vertices.end(); it++)
		{
			xmin = min(xmin, it->x);
			xmax = max(xmax, it->x);
			ymin = min(ymin, it->y);
			ymax = max(ymax, it->y);
		}

		// build a bounding triangle
		float halfWidth = xmax - xmin;
		float midX = (xmin + xmax) * 0.5f;
		float xminExt = midX - halfWidth * 1.38f,
			  xmaxExt = midX + halfWidth * 1.5f,
			  ymaxExt = ymin + 2.3f * (ymax - ymin);

		vertices.push_back(std::move(linalg::float2(xminExt, ymin - halfWidth * 0.12f)));
		auto i3 = --vertices.end();
		vertices.push_back(std::move(linalg::float2(xmaxExt, ymin - halfWidth * 0.18f)));
		auto i2 = --vertices.end();
		vertices.push_back(std::move(linalg::float2(midX, ymaxExt)));
		auto i1 = --vertices.end();
		
		edges.push_back(Edge(*i3, *i2));
		edges.push_back(Edge(*i2, *i1));
		edges.push_back(Edge(*i1, *i3));

		EdgeNode it = edges.begin();
		boundTriangle.pEdge[0] = it++;
		boundTriangle.pEdge[1] = it++;
		boundTriangle.pEdge[2] = it++;

		triangles.push_back(boundTriangle);
		it = edges.begin();
		for (int i = 0; i < 3; i++, it++)
		{
			it->pTriangle[0] = --triangles.end();
			it->pTriangle[1] = triangles.end();
		}
	}
	
	TriangleNode findContainingTriangle(const linalg::float2& v) // const
	{
		for (TriangleNode it = triangles.begin(); it != triangles.end(); it++)
		{
			if (pointInsideTriangle(v, *it)) 
				return it;
		}

		throw "no triangle contains this vertex";
	}

	void splitTriangle(TriangleNode triangleNode, const linalg::float2& v)
	{
		linalg::float2 verts[3];
		getTriangleVertices(*triangleNode, verts);

		EdgeNode e[3];
		edges.push_back(std::move(Edge(verts[0], v))); e[0] = --edges.end();
		edges.push_back(std::move(Edge(verts[1], v))); e[1] = --edges.end();
		edges.push_back(std::move(Edge(verts[2], v))); e[2] = --edges.end();

		for (int i = 0; i < 3; i++)
		{
			EdgeNode oldEdge = findEdge(*triangleNode, verts[i], verts[(i + 1) % 3]);
			Triangle triangle;
			triangle.pEdge[0] = e[i];
			triangle.pEdge[1] = e[(i + 1) % 3];
			triangle.pEdge[2] = oldEdge;
			triangles.push_back(triangle);
			e[i]->pTriangle[0] = --triangles.end();
			e[(i + 1) % 3]->pTriangle[1] = --triangles.end();

			if (oldEdge->pTriangle[0] == triangleNode)
				oldEdge->pTriangle[0] = --triangles.end();
			else
				oldEdge->pTriangle[1] = --triangles.end();
		}

		triangles.erase(triangleNode);
	}

	void flipEdge(EdgeNode edgeNode, EdgeNode affectedEdgeNodes[4])
	{
		TriangleNode t1 = edgeNode->pTriangle[0];
		linalg::float2 tv1[3];
		getTriangleVertices(*t1, tv1);
		
		TriangleNode t2 = edgeNode->pTriangle[1];
		linalg::float2 tv2[3];
		getTriangleVertices(*t2, tv2);

		linalg::float2 old_v1 = edgeNode->a;
		linalg::float2 old_v2 = edgeNode->b;
		linalg::float2 new_v1 = getOpposite(*edgeNode, *t1);
		linalg::float2 new_v2 = getOpposite(*edgeNode, *t2);

		Edge newEdge(new_v1, new_v2);
		edges.push_back(newEdge);
		EdgeNode newNode = --edges.end();
		
		EdgeNode e11 = findEdge(*t1, new_v1, old_v1);
		EdgeNode e12 = findEdge(*t2, new_v2, old_v1);
		EdgeNode e21 = findEdge(*t2, new_v2, old_v2);
		EdgeNode e22 = findEdge(*t1, new_v1, old_v2);
		t1->pEdge[0] = e11; t1->pEdge[1] = e12; t1->pEdge[2] = newNode;
		t2->pEdge[0] = e21; t2->pEdge[1] = e22; t2->pEdge[2] = newNode;
		newNode->pTriangle[0] = t1;
		newNode->pTriangle[1] = t2;
		
		if (e12->pTriangle[0] == t2) e12->pTriangle[0] = t1;
		else e12->pTriangle[1] = t1;

		if (e22->pTriangle[0] == t1) e22->pTriangle[0] = t2;
		else e22->pTriangle[1] = t2;


		edges.erase(edgeNode);

		affectedEdgeNodes[0] = e11;
		affectedEdgeNodes[1] = e12;
		affectedEdgeNodes[2] = e21;
		affectedEdgeNodes[3] = e22;
	}

	bool shouldFlip(EdgeNode edgeNode)
	{
		TriangleNode t1 = edgeNode->pTriangle[0];
		TriangleNode t2 = edgeNode->pTriangle[1];
		if (t2 == triangles.end()) // edge of boundary, can not flip
			return false;

		linalg::float2 o1 = getOpposite(*edgeNode, *t1);
		linalg::float2 v11 = edgeNode->a - o1;
		linalg::float2 v12 = edgeNode->b - o1;
		float cosA = linalg::dot(v11, v12) / (linalg::length(v11) * linalg::length(v12));
		float sinA = sqrtf(1.f - cosA * cosA);

		linalg::float2 o2 = getOpposite(*edgeNode, *t2);
		linalg::float2 v21 = edgeNode->a - o2;
		linalg::float2 v22 = edgeNode->b - o2;
		float cosB = linalg::dot(v21, v22) / (linalg::length(v21) * linalg::length(v22));
		float sinB = sqrtf(1.f - cosB * cosB);

		// sin(A+B) < 0  ==>  180 < A+B < 360
		return sinA * cosB + cosA * sinB < 0.f;
	}

	void addVertex(const linalg::float2& v)
	{
		std::queue<EdgeNode> queue;

		auto tNode = findContainingTriangle(v);
		queue.push(tNode->pEdge[0]);
		queue.push(tNode->pEdge[1]);
		queue.push(tNode->pEdge[2]);
		splitTriangle(tNode, v);

		while (!queue.empty())
		{
			auto eNode = queue.front(); queue.pop();
			if (shouldFlip(eNode))
			{
				EdgeNode affectedEdgeNodes[4];
				flipEdge(eNode, affectedEdgeNodes);
				for (const auto& en : affectedEdgeNodes)
					queue.push(en);
			}
		}
	}

	void getTriangleVertices(const Triangle& t, linalg::float2 tri[3]) const
	{
		tri[0] = t.pEdge[0]->a;
		tri[1] = t.pEdge[0]->b;
		tri[2] = t.pEdge[1]->a;
		if (tri[2] == tri[0] || tri[2] == tri[1])
			tri[2] = t.pEdge[1]->b;
	}

	bool pointInsideTriangle(const linalg::float2& v, const Triangle& t) const
	{
		linalg::float2 tri[3];
		getTriangleVertices(t, tri);
		float c1 = cross2(tri[1] - tri[0], v - tri[1]);
		float c2 = cross2(tri[2] - tri[1], v - tri[2]);
		float c3 = cross2(tri[0] - tri[2], v - tri[0]);
		return c1 * c2 > 0 && c1 * c3 > 0;
	}

	float cross2(const linalg::float2 a, const linalg::float2 b) const
	{
		return a.x * b.y - a.y * b.x;
	}

	linalg::float2 getOpposite(const Edge& e, const Triangle& t)
	{
		linalg::float2 verts[3];
		getTriangleVertices(t, verts);
		for (int i = 0; i < 3; i++)
		{
			if (verts[i] != e.a && verts[i] != e.b)
				return verts[i];
		}
		// error
		return linalg::float2(0, 0);
	}

	EdgeNode findEdge(const Triangle& t, const linalg::float2& a, const linalg::float2& b)
	{
		for (int i = 0; i < 3; i++)
		{
			if (t.pEdge[i]->a == a && t.pEdge[i]->b == b) return t.pEdge[i];
			if (t.pEdge[i]->a == b && t.pEdge[i]->b == a) return t.pEdge[i];
		}
		throw new std::exception("Exception while finding corresponding edge.");
		return edges.end();
	}


private:
	std::list<linalg::float2> vertices;
	std::list<Edge> edges;
	std::list<Triangle> triangles;
	Triangle boundTriangle;

	std::queue<Edge> edgeQueue;

};


int main()
{
	sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay Triangulation");
	Context context;
	context.pWindow = &window;

	ProblemConfig config;
	config.numPoints = 42;
	config.spaceWidth = 400;
	config.spaceHeight = 400;
	config.spaceX = 200;
	config.spaceY = 100;
	
	Problem problem(context, config);
	problem.start();
	DelaunayTriangulation dt;

	auto result = dt.Triangulate(problem.getVertices());
	if (result.error)
	{
		std::cerr << "Triangulation Error" << std::endl;
		exit(1);
	}

	dt.TriangulateSetup(problem.getVertices());

	std::vector<sf::Vertex> lines;
	for (const Edge& e : result.edges)
	{
		lines.push_back(sf::Vertex(sf::Vector2f(e.a.x, e.a.y)));
		lines.push_back(sf::Vertex(sf::Vector2f(e.b.x, e.b.y)));
	}

	bool isKeyPressedLastFrame = false;
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		//if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
		//{
		//	if (!isKeyPressedLastFrame)
		//	{
		//		auto result = dt.TriangulateStep();
		//		lines.clear();
		//		for (const Edge& e : result.edges)
		//		{
		//			lines.push_back(sf::Vertex(sf::Vector2f(e.a.x, e.a.y)));
		//			lines.push_back(sf::Vertex(sf::Vector2f(e.b.x, e.b.y)));
		//		}
		//	}
		//	isKeyPressedLastFrame = true;
		//}
		//else
		//{
		//	isKeyPressedLastFrame = false;
		//}

		window.clear(sf::Color::Black);
		problem.drawGeneratedVertices();
		if(lines.size() > 0)
			window.draw(&lines[0], lines.size(), sf::Lines);
		window.display();
	}

	return 0;
}