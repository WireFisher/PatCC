#ifndef PDLN_TRIANGLE_H
#define PDLN_TRIANGLE_H

#include "common_utils.h"
#ifdef OPENCV
#include "opencv2/opencv.hpp"
#endif


class Triangle;
class Triangle_inline;

class Point
{
    public:
        double x;
        double y;
        int    id;
        int    next;
        int    prev;

        Point();
        Point(double, double, int id=-1);
        ~Point();
        double calculate_distance(const Point*) const;
        double calculate_distance(double, double) const;
        int position_to_edge(const Point*, const Point*) const;
        int position_to_triangle(const Triangle*) const;
        int position_to_triangle(const Triangle_inline*) const;
        int is_in_region(double min_x, double max_x, double min_y, double max_y) const;
};


class Edge
{
    private:
        Point* head;
        Point* tail;                      /* the tail of this edge, constant */
        Edge*  twin_edge;                 /* the twin_edge edge, whose tail is the head of this edge and head is the tail of this edge */
        Edge*  next_edge_in_triangle;     /* the next_edge_in_triangle edge, whose tail is the head of this edge but head isn't the tail of this edge */
        Edge*  prev_edge_in_triangle;     /* the prev_edge_in_triangle edge, whose head is the tail of this edge but tail isn't the head of this edge */
        int    ref_count;
        Triangle* triangle;               /* the triangle which is composed by this edge and its next_edge_in_triangle and prev_edge_in_triangle */

    public:
        Edge(Point *head, Point *tail);
        ~Edge();
        Edge *generate_twins_edge();

        inline void ref_inc() {ref_count++;};
#ifdef OPENCV
        friend void draw_line(cv::Mat, Edge*, double, double, double, double, cv::Scalar);
#endif
        friend class Triangle;
        friend class Delaunay_Voronoi;
};


class Triangle
{
    private:
        Point* v[3];    /* vertexes of triangle */
        Edge*  edge[3];
        bool   is_leaf;
        bool   is_cyclic;
        int    remained_points_head;
        int    remained_points_tail;
        double circum_center[2];
        double circum_radius;
        bool   contain_virtual_polar;
        int    stack_ref_count;

        int  circum_circle_contains(Point*, double tolerance=PDLN_FLOAT_EQ_ERROR);
        bool really_on_circum_circle(Point *, double);

    public:
        Triangle();
        ~Triangle();
        double area();
        void get_center_coordinates();
        int find_best_candidate_point(Point*) const;
        bool contain_vertex(Point*);
        void calulate_circum_circle();

        int find_dividing_point(Point*);
        void set_remained_points(int, int);
        Point* pop_tail(Point*);

        friend class Delaunay_Voronoi;
        friend class Point;
        friend void  plot_triangles_into_file(const char *filename, std::vector<Triangle*>);
        friend class Triangle_pool;
};


class Triangle_inline
{
    public:
        Point v[3];
        bool is_cyclic;
        Triangle_inline() {};
        Triangle_inline(Point, Point, Point, bool = false);
        void check_cyclic();
        friend bool operator == (Triangle_inline, Triangle_inline);
};

#endif
