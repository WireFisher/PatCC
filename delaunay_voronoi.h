/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Mr. Yufeng Zhou,
  *  and then upgraded and merged into CoR by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef __DELAUNAY_VORONOI_H__
#define __DELAUNAY_VORONOI_H__

#include <vector>
#include <list>
#include <cmath>
#include <iostream>

#ifdef UNITTEST
#include "gtest/gtest_prod.h"
#include "opencv2/opencv.hpp"
#endif

using std::vector;

class Edge;
class Point;
class Triangle;


class Point
{
    private:
#ifdef UNITTEST
        FRIEND_TEST(DelaunayTriangulationTest, Random);
        friend void draw_line(cv::Mat, Edge*, double, double, double, double);
#endif
        friend class Delaunay_Voronoi;
        friend class Triangle;
        friend double det(const Point *, const Point *, const Point *);
        friend Point operator-(Point, Point);
        double x;    
        double y;
        int id;
        Triangle *current_triangle;

    public:
        Point();
        Point(double, double, int);
        ~Point();
        double calculate_distance(const Point*) const;
        int position_to_edge(const Point*, const Point*) const;
        int position_to_triangle(const Triangle*) const;
};


class Edge
{
    private:
#ifdef UNITTEST
        FRIEND_TEST(DelaunayTriangulationTest, Random);
        friend void draw_line(cv::Mat, Edge*, double, double, double, double);
#endif
        friend class Triangle;
        friend class Delaunay_Voronoi;
        Point *head;
        Point *tail;    /* the tail of this edge, constant */
        Edge *twin_edge;            /* the twin_edge edge, whose tail is the head of this edge and head is the tail of this edge */
        Edge *next_edge_in_triangle;            /* the next_edge_in_triangle edge, whose tail is the head of this edge but head isn't the tail of this edge */
        Edge *prev_edge_in_triangle;            /* the prev_edge_in_triangle edge, whose head is the tail of this edge but tail isn't the head of this edge */
        Triangle *triangle; /* the triangle which is composed by this edge and its next_edge_in_triangle and prev_edge_in_triangle */

    public:
        Edge(Point *head, Point *tail);
        ~Edge();
        Edge *generate_twins_edge();

};


class Triangle
{
    private:
        friend class Delaunay_Voronoi;
        friend class Point;
        Point *v[3];    /* vertexes of triangle */
        //Point center;    /* circumcenter */
        Edge *edge[3];
        bool is_leaf;
        int reference_count;    /* reference count, used to destruct */
        vector<Point*> remained_points_in_triangle;
        vector<Triangle*> children;
        double circum_center[2];
        double circum_radius;
        bool circum_circle_contains(Point*);

        Triangle(const Triangle &triangle);
        void initialize_triangle_with_edges(Edge*, Edge*, Edge*);
        void calulate_circum_circle();

    public:
        Triangle();
        Triangle(Point*, Point*, Point*);
        Triangle(Edge*, Edge*, Edge*);
        ~Triangle();
        void get_center_coordinates();
        int find_best_candidate_point();
        void check_and_set_twin_edge_relationship(Triangle*);
        bool contain_vertex(Point*);

};


struct Cell
{
    Point *center;
    vector<double> vertexes_lons;
    vector<double> vertexes_lats;
};


class Delaunay_Voronoi
{
    private:
        Cell *cells;
        vector<Triangle*> result_leaf_triangles;
        vector<Triangle*> triangle_pool;
        vector<Edge*> edge_pool;
        bool is_global_grid;
        int num_cells;
        Point *virtual_point[3];

        void check_and_set_twin_edge_relationship(vector<Triangle*>*);
        Point *generate_boundary_point(double, double, Triangle*, bool);
        void generate_initial_triangles(Triangle*, vector<Point*>*, vector<Point*>*, bool);
        void triangularization_process(Triangle*);
        void distribute_points_into_triangles(vector<Point*>*, vector<Triangle*>*);
        Triangle *search_triangle_with_point(Triangle*, const Point *pt);
        void generate_Voronoi_diagram();
        void extract_vertex_coordinate_values(int, bool, double**, double**, int*);
        void get_convex_set(int, double*, double*, double, double, int &, int **);
        Triangle* initialize_super_triangle(int, double*, double*, bool*);
        void clear_triangle_containing_virtual_point();

    public:
        Delaunay_Voronoi(int, double*, double*, bool, double, double, double, double, bool*);
        ~Delaunay_Voronoi();
        static bool is_triangle_legal(const Point *pt, const Edge *edge);
        void legalize_triangles(Point *pt, Edge *edge, vector<Triangle*>*);
        Edge *allocate_edge(Point *head, Point *tail);
        Triangle *allocate_Triangle(Point*, Point*, Point*);
        Triangle *allocate_Triangle(Edge*, Edge*, Edge*);
        vector<Edge*> get_all_delaunay_edge();
};

#endif
