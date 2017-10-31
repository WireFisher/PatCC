#ifndef H_DELAUNAY_2D
#define H_DELAUNAY_2D

#include <vector>

template <typename T>
class Vector2
{
	public:
		Vector2();
		Vector2(T _x, T _y);
		Vector2(const Vector2 &v);
		void set(const Vector2 &v);

		T dist2(const Vector2 &v);
		float dist(const Vector2 &v);

		T x;
		T y;
};


template <class T>
class Edge
{
	public:
		using VertexType = Vector2<T>;
		
		Edge(const VertexType &p1, const VertexType &p2) : p1(p1), p2(p2) {};
		Edge(const Edge &e) : p1(e.p1), p2(e.p2) {};

		VertexType p1;
		VertexType p2;
};


template <class T>
class Triangle
{
	public:
		using EdgeType = Edge<T>;
		using VertexType = Vector2<T>;
		
		Triangle(const VertexType &_p1, const VertexType &_p2, const VertexType &_p3):	p1(_p1), p2(_p2), p3(_p3), e1(_p1, _p2), 
                                                                                        e2(_p2, _p3), e3(_p3, _p1) {};
	
		bool containsVertex(const VertexType &v) { return p1 == v || p2 == v || p3 == v; };
		
		bool circumCircleContains(const VertexType &v);
	
		VertexType p1;
		VertexType p2;
		VertexType p3;
		EdgeType e1;				
		EdgeType e2;
		EdgeType e3;
};


template <class T>
class Delaunay
{
	public:
		using TriangleType = Triangle<T>;
		using EdgeType = Edge<T>;
		using VertexType = Vector2<T>;
		
		const std::vector<TriangleType>& triangulate(std::vector<VertexType> &vertices);
		
		const std::vector<TriangleType>& getTriangles() const { return _triangles; };
		const std::vector<EdgeType>& getEdges() const { return _edges; };
		const std::vector<VertexType>& getVertices() const { return _vertices; };

	private:
		std::vector<TriangleType> _triangles;
		std::vector<EdgeType> _edges;
		std::vector<VertexType> _vertices;
};

#endif
