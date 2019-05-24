/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


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
        int    id:31;
        int    mask:1;
        int    next;
        int    prev;

        Point();
        Point(double, double);
        Point(double, double, int, int = -1, int = -1);
        Point(double, double, int, bool, int = -1, int = -1);
        Point(PAT_REAL, PAT_REAL, int, bool, int = -1, int = -1);
        ~Point();
        double calculate_distance(const Point*) const;
        double calculate_distance(double, double) const;
        int position_to_edge(const Point*, const Point*) const;
        int position_to_triangle(const Point*, const Point*, const Point*) const;
        int position_to_triangle(const Triangle_inline*) const;
        int is_in_region(double min_x, double max_x, double min_y, double max_y) const;
};


class Edge
{
    private:
        int   head;
        int   tail;                      /* the tail of this edge, constant */
        Edge* twin_edge;                 /* the twin_edge edge, whose tail is the head of this edge and head is the tail of this edge */
        Edge* next_edge_in_triangle;     /* the next_edge_in_triangle edge, whose tail is the head of this edge but head isn't the tail of this edge */
        Edge* prev_edge_in_triangle;     /* the prev_edge_in_triangle edge, whose head is the tail of this edge but tail isn't the head of this edge */
        int   ref_count;
        Triangle* triangle;               /* the triangle which is composed by this edge and its next_edge_in_triangle and prev_edge_in_triangle */

    public:
        Edge();
        ~Edge();
        Edge *generate_twins_edge();

        inline void ref_inc() {ref_count++;};
#ifdef OPENCV
        friend void draw_line(cv::Mat, Point*, Edge*, double, double, double, double, cv::Scalar);
#endif
        friend class Triangle;
        friend class Delaunay_Voronoi;
};


class Triangle
{
    private:
        int      v[3];    /* index of vertexes */
        Edge*    edge[3];
        unsigned is_leaf:1;
        unsigned is_cyclic:1;
        unsigned is_virtual:1;
        int      remained_points_head;
        int      remained_points_tail;
        PAT_REAL circum_center[2];
        PAT_REAL circum_radius2;
        int      stack_ref_count;

        int  circum_circle_position_to(Point*, double=PDLN_ABS_TOLERANCE);

    public:
        Triangle();
        ~Triangle();
        void get_center_coordinates();
        int find_best_candidate_point(Point*);
        bool contain_vertex(int);
        void calulate_circum_circle(const Point*, const Point*, const Point*);

        int find_dividing_point(Point*);
        void set_remained_points(int, int);
        int pop_tail(Point*);

        friend class Delaunay_Voronoi;
        friend class Point;
        friend void  plot_triangles_into_file(const char *filename, std::vector<Triangle*>, Point*);
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
