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


template <typename T>
class Edge
{
	public:
		Edge(const Vector2<T> &p1, const Vector2<T> &p2) : p1(p1), p2(p2) {};
		Edge(const Edge &e) : p1(e.p1), p2(e.p2) {};

		Vector2<T> p1;
		Vector2<T> p2;
};

template <typename T>
class Triangle
{
	public:
		Triangle(const Vector2<T> &_p1, const Vector2<T> &_p2, const Vector2<T> &_p3):	p1(_p1), p2(_p2), p3(_p3), e1(_p1, _p2), 
                                                                                        e2(_p2, _p3), e3(_p3, _p1) {};
	
		bool containsVertex(const Vector2<T> &v) { return p1 == v || p2 == v || p3 == v; };
		
		bool circumCircleContains(const Vector2<T> &v);
	
		Vector2<T> p1;
		Vector2<T> p2;
		Vector2<T> p3;
		Edge<T> e1;				
		Edge<T> e2;
		Edge<T> e3;
};


template <typename T>
class Delaunay
{
	public:
		const std::vector<Triangle<T> >& triangulate(std::vector<Vector2<T> > &vertices);
		
		const std::vector<Triangle<T> >& getTriangles() const { return _triangles; };
		const std::vector<Edge<T> >& getEdges() const { return _edges; };
		const std::vector<Vector2<T> >& getVertices() const { return _vertices; };

	private:
		std::vector<Triangle<T> > _triangles;
		std::vector<Edge<T> > _edges;
		std::vector<Vector2<T> > _vertices;
};

#endif
