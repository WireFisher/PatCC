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
#include "dependency/MemoryPool/C-98/MemoryPool.h"

#ifdef UNITTEST
#include "gtest/gtest_prod.h"
#include "opencv2/opencv.hpp"
#endif

#define PDLN_CHECKSUM_FALSE (0xFFFFFFFF)

using std::vector;

class Edge;
class Point;
class Triangle;
class Triangle_Transport;

void sort_triangles(Triangle_Transport*, int);

class Point
{
    public:
#ifdef UNITTEST
        //friend void draw_line(cv::Mat, Edge*, double, double, double, double, cv::Scalar);
#endif
        //friend class Delaunay_Voronoi;
        //friend class Triangle;
        //friend class Delaunay_grid_decomposition;
        //friend double det(const Point *, const Point *, const Point *);
        //friend Point operator-(Point, Point);
        //friend bool operator == (Point, Point);
        double x;    
        double y;
        int id;
        Triangle *current_triangle;

    public:
        Point();
        Point(double, double, int id=-1);
        ~Point();
        double calculate_distance(const Point*) const;
        int position_to_edge(const Point*, const Point*) const;
        int position_to_triangle(const Triangle*) const;
        int is_in_region(double min_x, double max_x, double min_y, double max_y) const;
};


class Edge
{
    private:
#ifdef UNITTEST
        friend void draw_line(cv::Mat, Edge*, double, double, double, double, cv::Scalar);
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
        friend void plot_triangles_info_file(const char *filename, std::vector<Triangle*>);
        Point *v[3];    /* vertexes of triangle */
        //Point center;    /* circumcenter */
        Edge *edge[3];
        bool is_leaf;
        int reference_count;    /* reference count, used to destruct */
        vector<Point*> remained_points_in_triangle;
        vector<Triangle*> children;
        double circum_center[2];
        double circum_radius;
        int circum_circle_contains(Point*);

        void initialize_triangle_with_edges(Edge*, Edge*, Edge*);

    public:
        Triangle();
        Triangle(Edge*, Edge*, Edge*);
        ~Triangle();
        void get_center_coordinates();
        int find_best_candidate_point();
        void check_and_set_twin_edge_relationship(Triangle*);
        bool contain_vertex(Point*);
        void calulate_circum_circle();

};


struct Cell
{
    Point *center;
    vector<double> vertexes_lons;
    vector<double> vertexes_lats;
};

bool have_redundent_points(const double*, const double*, int);
void delete_redundent_points(double *&x, double *&y, int &num);

class Delaunay_Voronoi
{
    private:
        Cell *cells;
        vector<Triangle*> result_leaf_triangles;
        vector<Triangle*> triangle_pool;
        vector<Edge*> edge_pool;
        MemoryPool<Triangle, 0x10000*sizeof(Triangle), Edge*, Edge*, Edge*> triangle_allocator;
        MemoryPool<Edge, 0x10000*sizeof(Edge)> edge_allocator;
        bool is_global_grid;
        int num_cells;
        Point *virtual_point[3];
        int *global_index;

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
        bool is_angle_too_large(const Point *pt, const Edge *edge);
        bool is_angle_ambiguous(const Point *pt, const Edge *edge);
        const Point *get_lowest_point_of_four(const Point *, const Point *, const Point *, const Point *);
        double calculate_angle(const Point *, const Point *, const Point *);
        std::vector<Triangle*> find_triangles_intersecting_with_segment(Point, Point);

        bool is_triangle_legal(const Point *pt, const Edge *edge);
        bool is_triangle_legal(const Triangle *);
        bool is_triangle_ambiguous(const Point *pt, const Edge *edge);
        void relegalize_triangles(Point*, Edge*);
        void remove_leaf_triangle(Triangle*);
        unsigned calculate_triangles_checksum(Triangle_Transport*, int);

    public:
        Delaunay_Voronoi(int, double*, double*, int*, bool, double, double, double, double, bool*);
        ~Delaunay_Voronoi( );
        void legalize_triangles(Point *pt, Edge *edge, vector<Triangle*>*);
        Edge *allocate_edge(Point *head, Point *tail);
        Triangle *allocate_Triangle(Edge*, Edge*, Edge*);
        vector<Edge*> get_all_delaunay_edge();
        vector<Edge*> get_all_legal_delaunay_edge();
        bool is_all_leaf_triangle_legal();
        bool check_if_all_outer_edge_out_of_region(double, double, double, double);
        void get_triangles_in_region(double, double, double, double, Triangle_Transport *, int *, int);
        void update_all_points_coord(double *, double *, int);
        std::vector<Triangle*> search_cyclic_triangles_for_rotated_grid(Point, Point);
        void correct_cyclic_triangles(std::vector<Triangle*>, bool);
        void relegalize_all_triangles();
        void remove_triangles_on_or_out_of_boundary(double, double, double, double);
        unsigned calculate_triangles_intersected_checksum(Point, Point);   
        void get_triangles_intersecting_with_segment(Point, Point, Triangle_Transport*, int*, int);

        /* debug */
        void plot_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
};

class Triangle_Transport
{
    public:
        Point v[3];
        Triangle_Transport() {};
        Triangle_Transport(Point, Point, Point);
        friend bool operator == (Triangle_Transport, Triangle_Transport);
};

void plot_triangles_info_file(const char *filename, Triangle_Transport *t, int num);
void plot_triangles_info_file(const char *filename, std::vector<Triangle*>);
void save_triangles_info_file(const char *filename, Triangle_Transport *t, int num);

#endif
