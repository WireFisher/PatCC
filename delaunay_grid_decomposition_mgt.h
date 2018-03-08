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
    bool operator== (Boundry& boundry);
    bool operator!= (Boundry& boundry);
    bool operator<= (Boundry &boundry);
    void legalize();
    void legalize(const Boundry*, bool);
    void max(double, double, double, double);
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
    Boundry *real_boundry;
    double center[2];
    Midline midline;
    double *points_coord[2]; //(clean after having children)
    double *projected_coord[2];
    int *points_global_index;
    int len_expanded_points_coord_buf;
    int num_kernel_points;
    int num_expanded_points;
    int num_rotated_points; /* kernel + part of expanded */
    int node_type;
    bool non_monotonic;
    
    vector<int> processing_units_id;
    vector<pair<Search_tree_node*, bool> > neighbors;
    Delaunay_Voronoi *triangulation;
    void calculate_real_boundary();
    void fix_expanded_boundry(int index, int count);

public:    
    Search_tree_node(Search_tree_node*, double**, int*, int, Boundry, int type);
    ~Search_tree_node();
    void decompose_by_processing_units_number(double*, double**, int**, int*, Boundry*, vector<int>*, int);
    void decompose_by_fixed_longitude(double, double*, double**, int**, int*, Boundry*, vector<int>*);
    void split_local_points(Midline, double**, int**, int*);
    void split_processing_units_by_points_number(double*, int, int, vector<int>, vector<int>*);
    void update_processing_units_id(int);
    void update_processing_units_id(vector<int>);
    void generate_local_triangulation(bool);
    void generate_rotated_grid();
    void add_expanded_points(double **, int*, int);
    void add_neighbors(vector<Search_tree_node*>);
    bool check_if_all_outer_edge_out_of_kernel_boundry(Boundry *, bool);

    void search_points_in_halo(Boundry *inner_boundary, Boundry *outer_boundary, double *coord_values[2], int *global_idx, int *num_points_found);
    bool is_coordinate_in_halo(double x, double y, Boundry *inner, Boundry *outer);

    int get_num_kernel_points(){return num_kernel_points; };
    double** get_points_coord(){return points_coord; };
};

class Delaunay_grid_decomposition {
private:
    Processing_resource *processing_info;
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
    bool *active_processing_units_flag;
    //int* active_processing_common_id;
    double* workloads;
    
    void initialze_workload();
    int assign_polars(bool, bool);
    void decompose_common_node_recursively(Search_tree_node*);
    void decompose_with_fixed_longitude(double);
    void assign_cyclic_grid_for_single_processing_unit();
    bool have_local_processing_units_id(vector<int>);
    void update_workloads(int, vector<int>&);
    int expand_tree_node_boundry(Search_tree_node*, double);
    bool do_two_regions_overlap(Boundry, Boundry);
    Search_tree_node* alloc_search_tree_node(Search_tree_node*, double**, int*, int, Boundry, vector<int> &, int);

    void add_halo_points(Search_tree_node*, Boundry*, Boundry*);
    vector<Search_tree_node*> search_points_in_halo(Boundry*, Boundry*, double**, int*, int*);
    void search_down_for_points_in_halo(Search_tree_node*, Boundry*, Boundry*, vector<Search_tree_node*>&, double **, int*, int*);

    /* different decompositon consistency checking */
    bool check_leaf_node_triangulation_consistency(Search_tree_node*, int);
    void compute_common_boundry(Search_tree_node*, Search_tree_node*, Point*, Point*, Point*, Point*);
    
    /* process thread communication */
    int recv_triangles_from_remote(int, int, Triangle_Transport *, int, int);
    void send_triangles_to_remote(int, int, Triangle_Transport *, int, int);

    bool is_polar_node(Search_tree_node*) const;

    /* debug */
    void print_tree_node_info_recursively(Search_tree_node*);
    void save_ordered_triangles_into_file(Triangle_Transport *, int);

public:
    Delaunay_grid_decomposition(int, Processing_resource*, int);
    ~Delaunay_grid_decomposition();
    int generate_grid_decomposition();
    int generate_delaunay_triangulizition();
    vector<Search_tree_node*> get_local_leaf_nodes() {return local_leaf_nodes; };
    int generate_trianglulation_for_local_decomp();
    int generate_trianglulation_for_whole_grid();

    /* debug */
    void plot_local_triangles(const char*);
    void print_whole_search_tree_info();
    void merge_all_triangles();
};


class Triangle_ID_Only {
    public:
        int id[3];
        friend bool operator == (Triangle_ID_Only, Triangle_ID_Only);
};

class Grid_info_manager {
private:
    double *coord_values[2];
    int num_points;

    void gen_basic_grid();
    void gen_three_polar_grid();
    void gen_latlon_grid();
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
