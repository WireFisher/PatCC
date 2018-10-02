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
#include "common_utils.h"
#include "memory_pool.h"
#include "triangle.h"

#ifdef UNITTEST
#include "gtest/gtest_prod.h"
#endif

#define PDLN_CHECKSUM_FALSE (0xFFFFFFFF)

#define PDLN_UP     0
#define PDLN_LEFT   1
#define PDLN_DOWN   2
#define PDLN_RIGHT  3

using std::vector;

void sort_points_in_triangle(Triangle_pack*, int);
void sort_triangles(Triangle_pack*, int);

bool have_redundent_points(const double*, const double*, int);
void delete_redundent_points(double *&x, double *&y, int &num);

class Delaunay_Voronoi
{
    private:
        /* Storage */
        vector<Triangle*> all_leaf_triangles;
        Point*            all_points;

        /* Memory management */
        Triangle_pool triangle_allocator;
        Edge_pool edge_allocator;
        Triangle** triangle_stack;
        unsigned   stack_size;

        /* Property */
        bool   polar_mode;
        double tolerance;

        /* Grid info */
        int num_points;
        int vpolar_local_index;
        const double* x_ref;
        const double* y_ref;
        const int*    global_index;

        /* Triangulating stuff */
        Point* virtual_point[4];
        int*   point_idx_to_buf_idx;
        vector<Point*>    extra_virtual_point;
        vector<Triangle*> triangles_containing_vpolar;

        /* Consistency checking boundary */
        bool   have_bound;
        Point  bound_vertexes[4];
        double checking_threshold;
        vector<Triangle_pack> bound_triangles[4];

#ifdef DEBUG
        const double* x_store;
        const double* y_store;
#endif

        void triangulating_process(Triangle*, unsigned);
        void map_buffer_index_to_point_index();
        void push(unsigned *, Triangle*);

        void distribute_points_into_triangles(int, int, unsigned, unsigned);
        void link_remained_list(unsigned, unsigned, int*, int*);
        void swap_points(int, int);

        unsigned generate_initial_triangles();
        void clear_triangle_containing_virtual_point();
        bool is_angle_too_large(const Point *pt, const Edge *edge);
        bool is_angle_ambiguous(const Point *pt, const Edge *edge);
        const Point *get_lowest_point_of_four(const Point *, const Point *, const Point *, const Point *);
        Edge* generate_twins_edge(Edge*);

        bool is_triangle_legal(const Point *pt, const Edge *edge);
        bool is_triangle_legal(const Triangle *);
        bool is_triangle_ambiguous(const Point *pt, const Edge *edge);
        void relegalize_triangles(Point*, Edge*);
        void remove_leaf_triangle(Triangle*);
        void update_virtual_polar_info();

        Triangle_pack pack_triangle(Triangle*);
        void add_to_bound_triangles(Triangle_pack&, unsigned);

        void legalize_triangles(Point *pt, Edge *edge, unsigned, unsigned*);

        Edge* allocate_edge(Point *head, Point *tail);
        Triangle* allocate_Triangle(Edge*, Edge*, Edge*);
        void initialize_triangle_with_edges(Triangle*, Edge*, Edge*, Edge*, bool force=false);

        inline void ref_inc(Edge* e) {e->ref_count++;};
        inline void ref_dec(Edge* e) {
            e->ref_count--;
            if(e->ref_count <= 0) {
                if (e->twin_edge)
                    e->twin_edge->twin_edge = NULL;
                edge_allocator.deleteElement(e);
            }
        };
        void clean_triangle(Triangle*);

    public:
        Delaunay_Voronoi();
        ~Delaunay_Voronoi();

        void add_points(const double* x, const double* y, const int* global_idx, int num);
        void set_virtual_polar_index(int idx);
        void set_origin_coord(const double *x_origin, const double *y_origin);
        void triangulate();
        void set_checksum_bound(double, double, double, double, double);
        void set_polar_mode(bool);

        bool is_all_leaf_triangle_legal();
        vector<Edge*> get_all_delaunay_edge();
        vector<Edge*> get_all_legal_delaunay_edge();
        void get_triangles_in_region(double, double, double, double, Triangle_pack *, int *, int);
        void update_all_points_coord(double *, double *, int);
        void relegalize_all_triangles();
        void remove_triangles_on_or_out_of_boundary(double, double, double, double);

        unsigned calculate_checksum(Point, Point, double = 0);
        int bound_direction(const Point*, const Point*);
        unsigned cal_checksum(Point, Point, double = 0);

        bool is_triangle_in_circle(Triangle*, Point, double);

        void remove_triangles_in_circle(Point, double);
        void remove_triangles_on_segment(Point, Point);
        void recognize_cyclic_triangles();

        void update_points_coord_y(double, vector<int> *);
        void remove_triangles_only_containing_virtual_polar();
        void uncyclic_all_points();
        void remove_triangles_till(int);

        void make_final_triangle();
        void make_final_triangle_pack();
        void make_bounding_triangle_pack();

        /* debug */
        void save_original_points_into_file();
#ifdef OPENCV
        void plot_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_projection_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_original_points_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_current_step_into_file(const char*);
#endif
};

#ifdef OPENCV
void plot_triangles_into_file(const char *filename, Triangle_pack *t, int num, bool plot_cyclic_triangles=true);
void plot_triangles_into_file(const char *filename, std::vector<Triangle*>);
#endif
void save_triangles_info_file(const char *filename, Triangle_pack *t, int num);

#endif
