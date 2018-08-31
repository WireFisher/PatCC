/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#ifndef DELAUNAY_GRID_DECOMPOSITION_MGT_H
#define DELAUNAY_GRID_DECOMPOSITION_MGT_H

#include "processing_unit_mgt.h"
#include "delaunay_voronoi_2D.h"

#define PDLN_LON 0
#define PDLN_LAT 1

#define PDLN_FLOAT_EQ_ERROR (1e-10)

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
    Boundry& operator* (double);
    bool operator== (const Boundry& boundry) const;
    bool operator!= (const Boundry& boundry) const;
    bool operator<= (const Boundry &boundry) const;
    void legalize();
    void legalize(const Boundry*, bool);
    void max(double, double, double, double);
    void max(const Boundry);
    void squeeze(const Boundry*, double);
    void move_close(double**, int, int);
};

class Search_tree_node;
typedef vector<pair<Search_tree_node*, bool> > Neighbors;

class Search_tree_node {
private:
    Search_tree_node *parent;
    Search_tree_node *children[3];
    int node_type;
    int region_id;

    Boundry* kernel_boundry;
    Boundry* expanded_boundry;
    Boundry* rotated_kernel_boundry;
    Boundry* rotated_expanded_boundry;
    Boundry* real_boundry;

    double  center[2];
    double* points_coord[2]; //(clean after having children)
    double* projected_coord[2];
    int*    points_global_index;

    int     len_expanded_points_coord_buf;
    int     num_kernel_points;
    int     num_expanded_points;
    int     num_projected_points; /* number of points already been projected, kernel + part of expanded */
    bool    non_monotonic;
    Midline midline;

    vector<int> region_ids;
    int*        group_intervals;
    int         num_groups;

    int expanding_scale[4];
    Neighbors neighbors;
    Delaunay_Voronoi* triangulation;

    vector<int>* polars_local_index;
    int          virtual_point_local_index;
    
    /* consistency check */
    int num_neighbors_on_boundry[4];
    int edge_expanding_count[4];

    friend class Delaunay_grid_decomposition;

    void fix_view_point();
    void calculate_real_boundary();
    void fix_expanded_boundry(int index, int count);
    double load_polars_info();
    void reset_polars();

public:    
    Search_tree_node(Search_tree_node*, double**, int*, int, Boundry, int type);
    ~Search_tree_node();
    void decompose_by_processing_units_number(double*, double**, int**, int*, Boundry*, vector<int>*, int, int** =NULL, int* =NULL);
    void decompose_by_fixed_longitude(double, double*, double**, int**, int*, Boundry*, vector<int>*);
    void split_local_points(Midline, double**, int**, int*);
    void split_processing_units_by_points_number(double*, int, int, vector<int>, vector<int>*);
    void update_region_ids(int);
    void update_region_ids(vector<int>);
    void generate_local_triangulation(bool, int);
    void project_grid();
    void add_expanded_points(double *, double *, int*, int);
    void add_expanded_points(double **, int*, int);
    void add_neighbors(vector<Search_tree_node*>);
    void init_num_neighbors_on_boundry(int);
    bool expanding_success(Boundry *, Boundry *, bool);

    void reduce_num_neighbors_on_boundry(unsigned);
    void clear_expanding_count(unsigned);

    Boundry expand();
    void set_groups(int *, int);

    static void search_points_in_halo(const Boundry*, const Boundry*, double*const *, const int*, int, double**, int*, int*);
    void search_points_in_halo(const Boundry*, const Boundry*, double**, int*, int*);
    static bool is_coordinate_in_halo(double x, double y, const Boundry *inner, const Boundry *outer);
    static void count_points(double**, int*, int, int, Midline, int*);

    int get_num_kernel_points(){return num_kernel_points; };
    double** get_points_coord(){return points_coord; };
};

class Delaunay_grid_decomposition {
private:
    Search_tree_node*         search_tree_root;
    Search_tree_node*         current_tree_node;
    vector<Search_tree_node*> all_leaf_nodes;
    vector<Search_tree_node*> local_leaf_nodes;
    int min_points_per_chunk;

    /* grid info */
    int  original_grid;
    bool is_cyclic;
    int  num_points;
    int  num_inserted;

    /* proc info */
    Processing_resource* processing_info;
    bool*     active_processing_units_flag;
    double*   workloads;
    double    average_workload;
    int*      regionID_to_unitID;
    unsigned* regionID_to_checksum;
    int*      all_group_intervals;

    void initialze_workload();
    int assign_polars(bool, bool);
    void decompose_common_node_recursively(Search_tree_node*, bool =true);
    void decompose_with_fixed_longitude(double);
    void assign_cyclic_grid_for_single_processing_unit();
    bool have_local_region_ids(vector<int>);
    void update_workloads(int, vector<int>&);
    int expand_tree_node_boundry(Search_tree_node*, double);
    bool do_two_regions_overlap(Boundry, Boundry);
    Search_tree_node* alloc_search_tree_node(Search_tree_node*, double**, int*, int, Boundry, vector<int> &, int);
    int insert_virtual_points(double *coord_values[2], Boundry *boundry, int num_points);
    vector<Search_tree_node*> adjust_expanding_boundry(const Boundry*, Boundry*, double, double**, int*, int*);
    //void search_halo_points();
    static double adjust_subrectangle(double, double, double**, int*, int, int, Boundry*, int, int);

    void add_halo_points(Search_tree_node*, Boundry*, Boundry*);
    vector<Search_tree_node*> search_halo_points_from_top(const Boundry*, const Boundry*, double**, int*, int*);
    void search_down_for_points_in_halo(Search_tree_node*, const Boundry*, const Boundry*, vector<Search_tree_node*>&, double **, int*, int*);
    void search_halo_points_from_buf(Boundry*, Boundry*, double**, int*, int*);

    static void halo_to_rectangles(Boundry, Boundry, Boundry*);
    static void rectangles_to_halo(Boundry*, Boundry*);
    static int classify_points(double**, int*, int, Boundry, int);

    /* different decompositon consistency checking */
    bool check_leaf_node_triangulation_consistency(Search_tree_node*, int);
    unsigned compute_common_boundry(Search_tree_node*, Search_tree_node*, Point*, Point*, Point*, Point*);
    
    /* process thread communication */
    int recv_triangles_from_remote(int, int, Triangle_Transport *, int, int);
    void send_triangles_to_remote(int, int, Triangle_Transport *, int, int);

    bool is_polar_node(Search_tree_node*) const;

    /* debug */
    void print_tree_node_info_recursively(Search_tree_node*);
    void save_ordered_triangles_into_file(Triangle_Transport *&, int);

public:
    Delaunay_grid_decomposition(int, Processing_resource*, int);
    ~Delaunay_grid_decomposition();
    int generate_grid_decomposition(bool =true);
    int grid_full_decomposition();
    int generate_delaunay_triangulizition();
    vector<Search_tree_node*> get_local_leaf_nodes() {return local_leaf_nodes; };
    int generate_trianglulation_for_local_decomp();
    int generate_trianglulation_for_whole_grid();
    void send_recv_checksums_with_neighbors(Search_tree_node*, unsigned*, unsigned*, vector<MPI_Request*> *, int);
    bool are_checksums_identical(Search_tree_node*, unsigned*, unsigned*);

    /* debug */
    void print_whole_search_tree_info();
    void merge_all_triangles();
#ifdef NETCDF
    void plot_grid_decomposition(const char*);
    void plot_local_triangles(const char*);
#endif
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
    double min_lon;
    double max_lon;
    double min_lat;
    double max_lat;

    void gen_basic_grid();
#ifdef NETCDF
    void gen_three_polar_grid();
    void gen_latlon_grid();
    void gen_latlon_90_grid();
#endif

public:
    /* for unittest */
    Grid_info_manager();
    Grid_info_manager(double *coord[2], int num);
    virtual ~Grid_info_manager();
    virtual double** get_grid_coord_values(int);
    virtual int get_grid_num_points(int);
    virtual void get_grid_boundry(int, double*, double*, double*, double*);
    virtual void set_grid_boundry(int, double, double, double, double);
    virtual bool is_grid_cyclic(int);
};

extern Grid_info_manager *grid_info_mgr;
#endif
