/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#ifndef DELAUNAY_GRID_DECOMPOSITION_MGT_H
#define DELAUNAY_GRID_DECOMPOSITION_MGT_H

#define PDLN_LON 0
#define PDLN_LAT 1
#include "processing_unit_mgt.h"
#include "delaunay_voronoi_2D.h"

using std::pair;

struct Midline {
    int type;
    double value;
};

class Boundry {
public:
    double min_lon;
    double max_lon;
    double min_lat;
    double max_lat;

    Boundry() {};
    Boundry(double l, double r, double b, double t): min_lon(l), max_lon(r), min_lat(b), max_lat(t) {};
    Boundry& operator= (Boundry& boundry);
    Boundry& operator* (double);
    void legalize();
    void legalize(const Boundry*, bool);
};


class Search_tree_node {
private:
    friend class Delaunay_grid_decomposition;
    Search_tree_node *parent;
    Search_tree_node *children[3];
    Boundry *kernel_boundry;
    Boundry *expanded_boundry;
    Boundry *rotated_kernel_boundry;
    Boundry *rotated_expanded_boundry;
    double center[2];
    //double expanding_ratio;
    Midline midline;
    double *local_cells_coord[2]; //(clean after having children) //TODO:change name
    double *rotated_cells_coord[2];
    int *local_cells_global_index;
    int len_expanded_cells_coord_buf;
    int num_local_kernel_cells;//TODO:change name
    int num_local_expanded_cells;//TODO:change name
    int num_rotated_cells; /* kernel + part of expanded */
    int node_type;
    
    vector<int> processing_units_id;
    vector<pair<Search_tree_node*, bool> > neighbors;
    Delaunay_Voronoi *triangulation;

public:    
    Search_tree_node(Search_tree_node*, double**, int*, int, Boundry, int type=0); //FIXME: remove default value;
    ~Search_tree_node();
    void decompose_iteratively(double*, double **, int**, int*, Boundry*, vector<int>*, int);
    void decompose_with_certain_line(Midline, double**, int**, int*);
    void update_processing_units_id(int);
    void update_processing_units_id(vector<int>);
    void generate_local_triangulation(bool);
    void generate_rotated_grid();
    void add_expanded_points(double **, int*, int);
    void add_neighbors(vector<Search_tree_node*>);
    //bool check_expanded_triangle_consistency();
    bool check_if_all_outer_edge_out_of_kernel_boundry(Boundry *, bool);

    int get_num_local_kernel_cells(){return this->num_local_kernel_cells; };
    double** get_local_cells_coord(){return this->local_cells_coord; };
};

class Delaunay_grid_decomposition {
private:
    Processing_resource *processing_info;
    //Workload_info *workload_info;
    //Remap_grid_class *original_grid;
    int original_grid;
    Search_tree_node *search_tree_root;
    Search_tree_node *current_tree_node;
    vector<Search_tree_node*> local_leaf_nodes;
    //vector<delaunay_triangulation*> local_triangulations;
    //int *num_local_nodes_per_thread;
    //int num_total_nodes;
    //int *local_units_id;
    int min_num_points_per_chunk;
    bool is_cyclic;
    int* active_processing_common_id;
    double* workloads;
    
    int rotate_grid();
    void initialze_workload();
    int assign_polars(bool, bool);
    void decompose_common_node_recursively(Search_tree_node*);
    void assign_cyclic_grid_for_single_unit();
    bool have_local_processing_units_id(vector<int>);
    void update_workloads(int, vector<int>&);
    int expand_tree_node_boundry(Search_tree_node*, double);
    vector<Search_tree_node*> search_points_in_region(Boundry, double**, int*, int*);
    void transform_into_rectangle(Boundry, Boundry, Boundry*);
    void search_leaf_nodes_overlapping_with_region_recursively(Search_tree_node*, Boundry, vector<Search_tree_node*>&);
    bool do_two_regions_overlap(Boundry, Boundry);

    /* different decompositon consistency checking */
    bool check_leaf_node_triangulation_consistency(Search_tree_node*, int);
    void compute_common_boundry(Search_tree_node*, Search_tree_node*, Point*, Point*, Point*, Point*);
    bool check_triangles_consistency(Triangle_Transport*, Triangle_Transport*, int);
    //void get_triangles_intersecting_with_segment_from_remote(int, Point, Point, vector<Triangle_Transport>*);
    
    /* process thread communication */
    int recv_triangles_from_remote(int, int, Triangle_Transport *, int, int);
    void send_triangles_to_remote(int, int, Triangle_Transport *, int, int);

    /* debug */
    void print_tree_node_info_recursively(Search_tree_node*);
    void save_ordered_triangles_into_file(Triangle_Transport *, int);

public:
    //Delaunay_grid_decomposition(int, Processing_resource*);
    Delaunay_grid_decomposition(int, Processing_resource*, int);
    ~Delaunay_grid_decomposition();
    int generate_grid_decomposition();
    int generate_delaunay_triangulizition();
    vector<Search_tree_node*> get_local_leaf_nodes() {return this->local_leaf_nodes; };
    int generate_trianglulation_for_local_decomp();
    int generate_trianglulation_for_whole_grid();

    /* debug */
    void plot_local_triangles(const char*);
    void print_whole_search_tree_info();
    void merge_all_triangles();
};


class Grid_info_manager {
private:
    double *coord_values[2];
    int num_points;
public:
    /* for unittest */
    Grid_info_manager();
    Grid_info_manager(double *coord[2], int num);
    virtual ~Grid_info_manager();
    virtual double** get_grid_coord_values(int);
    virtual int get_grid_num_points(int);
    virtual void get_grid_boundry(int, double*, double*, double*, double*);
    virtual int get_polar_points(int, char);
    virtual bool is_grid_cyclic(int);
};

extern Grid_info_manager *grid_info_mgr;
#endif
