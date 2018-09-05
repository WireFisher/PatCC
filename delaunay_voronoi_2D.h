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

#include "processing_unit_mgt.h" // only for output debugging message
#include <vector>
#include <list>
#include <cmath>
#include <iostream>
#include "dependency/MemoryPool/C-98/MemoryPool.h"

#ifdef UNITTEST
#include "gtest/gtest_prod.h"
#endif

#ifdef OPENCV
#include "opencv2/opencv.hpp"
#endif

#define PDLN_CHECKSUM_FALSE (0xFFFFFFFF)

#define FLOAT_ERROR_LOW ((double) 1e-8)
#define FLOAT_ERROR ((double) 1e-10) // if less than 1e-10, will three point in a line, if more than 1e-15, will not pass check
#define FLOAT_ERROR_HI ((double) 1e-10) // normal grid less than 1e-11

#define PI ((double) 3.1415926535897932384626433)

using std::vector;

class Edge;
class Point;
class Triangle;
class Triangle_Transport;

void sort_points_in_triangle(Triangle_Transport*, int);
void sort_triangles(Triangle_Transport*, int);

class Point
{
    public:
        double x;    
        double y;
        int id;
        Triangle *current_triangle;

    public:
        Point();
        Point(double, double, int id=-1);
        ~Point();
        double calculate_distance(const Point*) const;
        double calculate_distance(double, double) const;
        int position_to_edge(const Point*, const Point*) const;
        int position_to_triangle(const Triangle*) const;
        int is_in_region(double min_x, double max_x, double min_y, double max_y) const;
};


class Edge
{
    private:
#ifdef OPENCV 
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
        friend void plot_triangles_into_file(const char *filename, std::vector<Triangle*>);
        Point *v[3];    /* vertexes of triangle */
        //Point center;    /* circumcenter */
        Edge *edge[3];
        bool is_leaf;
        bool is_cyclic;
        vector<Point*> remained_points_in_triangle;
        double circum_center[2];
        double circum_radius;
        int circum_circle_contains(Point*, double tolerance=FLOAT_ERROR);
        bool really_on_circum_circle(Point *, double);
        bool contain_virtual_polar;

        void initialize_triangle_with_edges(Edge*, Edge*, Edge*, bool force=false);

    public:
        Triangle();
        Triangle(Point*, Point*, Point*);
        Triangle(Edge*, Edge*, Edge*);
        ~Triangle();
        double area();
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
        Point *virtual_point[4];
        vector<Point*> extra_virtual_point;
        int *global_index;
        vector<Triangle*> triangles_containing_vpolar;
        double lat_nearest_vpolar;
        int vpolar_local_index;
        double tolerance;
        const double *x_ref;
        const double *y_ref;

        void check_and_set_twin_edge_relationship(vector<Triangle*>*);
        void triangularization_process(Triangle*);
        void distribute_points_into_triangles(vector<Point*>*, vector<Triangle*>*);
        vector<Triangle*> generate_initial_triangles(int, double*, double*, bool*);
        void clear_triangle_containing_virtual_point();
        bool is_angle_too_large(const Point *pt, const Edge *edge);
        bool is_angle_ambiguous(const Point *pt, const Edge *edge);
        const Point *get_lowest_point_of_four(const Point *, const Point *, const Point *, const Point *);
        double calculate_angle(const Point *, const Point *, const Point *);
        std::vector<Triangle*> find_triangles_intersecting_with_segment(Point, Point, double);
        Edge* generate_twins_edge(Edge*);

        bool is_triangle_legal(const Point *pt, const Edge *edge);
        bool is_triangle_legal(const Triangle *);
        bool is_triangle_ambiguous(const Point *pt, const Edge *edge);
        void relegalize_triangles(Point*, Edge*);
        void remove_leaf_triangle(Triangle*);
        unsigned calculate_triangles_checksum(Triangle_Transport*, int);
        void update_virtual_polar_info();

    public:
        Delaunay_Voronoi(int, double*, double*, const double*, const double*, int*, bool, double, double, double, double, bool*, int virtual_polar_local_index=-1);
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
        unsigned calculate_triangles_intersected_checksum(Point, Point, double threshold = 0);   
        void get_triangles_intersecting_with_segment(Point, Point, Triangle_Transport*, int*, int, double threshold = 0);
        bool is_triangle_in_circle(Triangle*, Point, double);
        void remove_triangles_in_circle(Point, double);
        void remove_triangles_on_segment(Point, Point);
        void recognize_cyclic_triangles();

        void update_points_coord_y(double, vector<int> *);
        void remove_triangles_only_containing_virtual_polar();
        void uncyclic_all_points();
        void remove_triangles_till(int);

        /* debug */
        void save_original_points_into_file();
#ifdef OPENCV
        void plot_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_projection_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_original_points_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_current_step_into_file(const char*);
#endif
};

class Triangle_Transport
{
    public:
        Point v[3];
        Triangle_Transport() {};
        Triangle_Transport(Point, Point, Point);
        friend bool operator == (Triangle_Transport, Triangle_Transport);
};

#ifdef OPENCV
void plot_triangles_into_file(const char *filename, Triangle_Transport *t, int num, bool plot_cyclic_triangles=true);
void plot_triangles_into_file(const char *filename, std::vector<Triangle*>);
#endif
void save_triangles_info_file(const char *filename, Triangle_Transport *t, int num);

#endif
