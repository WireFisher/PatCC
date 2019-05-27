/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DELAUNAY_GRID_DECOMPOSITION_MGT_H
#define DELAUNAY_GRID_DECOMPOSITION_MGT_H

#include "processing_unit_mgt.h"
#include "triangulation.h"

#define PDLN_LON 0
#define PDLN_LAT 1
#define PDLN_HIGH_BOUNDRY_SHIFTING (1e-6)


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
    void move_close(double**, int, int);
};


struct Grid_info
{
    double* coord_values[2];
    bool*   mask;
    int     num_total_points;
    int     num_vitual_poles;
    int     num_fence_points;
    Boundry boundary;
    bool    is_cyclic;
};

class Search_tree_node;
class Delaunay_grid_decomposition;
typedef vector<pair<Search_tree_node*, bool> > Neighbors;


class Search_tree_node {
private:
    Search_tree_node *parent;
    Search_tree_node *children[3];
    int  node_type;
    int  region_id;
    bool fast_triangulate;

    Boundry* kernel_boundry;
    Boundry* expand_boundry;
    Boundry* project_boundry;
    Boundry* real_boundry;

    double  center[2];
    double* kernel_coord[2];
    double* expand_coord[2];
    double* projected_coord[2];

    int*    kernel_index;
    int*    expand_index;
    bool*   kernel_mask;
    bool*   expand_mask;

    int     len_expand_coord_buf;
    int     num_kernel_points;
    int     num_expand_points;
    int     num_projected_points; /* number of points already been projected, kernel + part of expanded */
    int     num_old_points;
    bool    non_monotonic;
    Midline midline;

    int  ids_start;
    int  ids_end;
    int* group_intervals;
    int  num_groups;

    int expanding_scale[4];
    Neighbors neighbors;
    Delaunay_Voronoi* triangulation;
    int  bind_with;
    bool is_bind;

    vector<int>* polars_local_index;
    double       shifted_polar_lat;
    int          virtual_point_local_index;
    
    /* consistency check */
    int num_neighbors_on_boundry[4];
    int edge_expanding_count[4];

    void sort_by_line(Midline*, int*, int*);
    static void sort_by_line_internal(double**, int*, bool*, Midline*, int, int, int*, int*);

    void fix_view_point();
    void calculate_real_boundary();
    void fix_expand_boundry(int index, int count);
    double load_polars_info();
    void reset_polars(double*);
    void calculate_latitude_circle_projection(double, Point*, double*);
    void calculate_cyclic_boundary_projection(unsigned, Point*, Point*);

public:    
    Search_tree_node(Search_tree_node*, double**, int*, bool*, int, Boundry, int type);
    ~Search_tree_node();

    /* Grid Decomposition */
    void decompose_by_processing_units_number(double*, double**, int**, bool**, int*, Boundry*, int*, int*, int, int**, int*, int);
    void divide_at_fix_line(Midline, double**, int**, bool**, int*);
    void reorganize_kernel_points(double, double, double, double, int, int, Midline*, int*, int);
    static int divide_points(double**, int*, bool*, double, double, double, double, int, int, int, Midline*, int*, double, double, double, int);

    /* Getter & Setter */
    void update_region_ids(int, int);
    void set_groups(int *, int);
    inline int ids_size() {return ids_end - ids_start; };

    /* Triangulation */
    void project_grid();
    void generate_local_triangulation(bool, int, int, bool);

    /* Expanding */
    Boundry expand();
    void add_expand_points(double*, double*, int*, bool*, int);
    void add_expand_points(double**, int*, bool*, int);
    void add_neighbors(vector<Search_tree_node*>);
    void init_num_neighbors_on_boundry(int);
    bool expanding_success(bool*);

    /* Points searching */
    static void search_points_in_halo_internal(const Boundry*, const Boundry*, double *const *, const int*, const bool*, int, double**, int*, bool*, int*);
    void search_points_in_halo(const Boundry*, const Boundry*, double**, int*, bool*, int*);
    static bool is_coordinate_in_halo(double x, double y, const Boundry *inner, const Boundry *outer);

    /* Consistency checking */
    void reduce_num_neighbors_on_boundry(unsigned long);
    void clear_expanding_count(unsigned long);

    friend class Delaunay_grid_decomposition;
    friend bool node_ptr_comp(Search_tree_node*, Search_tree_node*);
    friend void decompose_common_node_recursively(Delaunay_grid_decomposition *, Search_tree_node *, int, bool);
    friend void extend_search_tree(Delaunay_grid_decomposition *, Search_tree_node *, const Boundry*, int, int);
};

class Delaunay_grid_decomposition {
public:
    Delaunay_grid_decomposition(Grid_info, Processing_resource*, int);
    ~Delaunay_grid_decomposition();

    int generate_grid_decomposition(bool =true);
    int generate_trianglulation_for_local_decomp();
    vector<Search_tree_node*> get_local_leaf_nodes() {return local_leaf_nodes; };

    /* Debug */
    void print_whole_search_tree_info();
    void merge_all_triangles(bool);

#ifdef OPENCV
    void plot_grid_decomposition(const char*);
    void plot_local_triangles(const char*);
#endif

private:
    /* Main processes */
    int  initialze_workloads(bool, bool);
    void initialze_buffer();
    int  assign_polars(bool, bool);

    /* Pre-treatment */
    int calculate_num_inserted_points(Boundry*, int);
    int dup_inserted_points(double**, bool**, Boundry*, int);

    /* Helper */
    bool have_local_region_ids(int, int);
    void update_workloads(int, int, int, bool);
    Search_tree_node* alloc_search_tree_node(Search_tree_node*, double**, int*, bool*, int, Boundry, int, int, int, bool=false);
    bool is_polar_node(Search_tree_node*) const;
    void set_binding_relationship();

    /* Grid Expanding */
    int expand_tree_node_boundry(Search_tree_node*, double);
    vector<Search_tree_node*> adjust_expanding_boundry(const Boundry*, Boundry*, double*, double**, int*, bool*, bool*, int*);
    bool do_two_regions_overlap(Boundry, Boundry);
    static void adjust_subrectangle(double, double, double**, int*, bool*, int, int, Boundry*, int, int, int*, int*);
    static int move_together(double**, int*, bool*, int*, int*, Boundry);
    static void halo_to_rectangles(Boundry, Boundry, Boundry*);
    static void rectangles_to_halo(Boundry*, Boundry*);
    static int classify_points(double**, int*, bool*, int, Boundry, int);

    /* Points searching */
    void search_down_for_points_in_halo(Search_tree_node*, const Boundry*, const Boundry*, vector<Search_tree_node*>*, double **, int*, bool*, int*);

    /* Consistency checking */
    bool check_leaf_node_triangulation_consistency(Search_tree_node*, int);
    unsigned compute_common_boundry(Search_tree_node*, Search_tree_node*, Point*, Point*, Point*, Point*);
    void send_recv_checksums_with_neighbors(Search_tree_node*, unsigned long*, unsigned long*, vector<MPI_Request*> *, int);
    bool are_checksums_identical(Search_tree_node*, unsigned long*, unsigned long*);
    void send_checksum_to_remote(int src_common_id, int dst_common_id, unsigned long* , int tag, MPI_Request** req);
    void recv_checksum_from_remote(int src_common_id, int dst_common_id, unsigned long*, int tag, MPI_Request** req);
    
    /* Process thread communication */
    int recv_triangles_from_remote(int, int, Triangle_inline *, int, int);
    void send_triangles_to_remote(int, int, Triangle_inline *, int, int);

    /* Debug */
    void print_tree_node_info_recursively(Search_tree_node*);
    void save_unique_triangles_into_file(Triangle_inline *&, int, bool);

    /* Search tree info */
    Search_tree_node*         search_tree_root;
    Search_tree_node*         current_tree_node;
    vector<Search_tree_node*> all_leaf_nodes;
    vector<Search_tree_node*> local_leaf_nodes;
    int  min_points_per_chunk;

    /* Grid info */
    int     original_grid;
    bool    is_cyclic;
    double* coord_values[2];
    bool*   mask;
    int*    global_index;
    int     num_points;
    int     num_fence_points;
    Boundry boundary_from_user;

    /* Proc info */
    Processing_resource* processing_info;
    int       num_regions;
    bool*     active_processing_units_flag;
    bool      is_local_proc_active;
    double*   workloads;
    double    average_workload;
    int*      regionID_to_unitID;
    int*      all_group_intervals;

    /* Temp buffer */
    double** buf_double[2];
    int**    buf_int;
    bool**    buf_bool;

    friend void decompose_common_node_recursively(Delaunay_grid_decomposition *, Search_tree_node *, int, bool);
    friend void extend_search_tree(Delaunay_grid_decomposition *, Search_tree_node *, const Boundry*, int, int);
};


class Triangle_ID_Only {
    public:
        int id[3];
        friend bool operator == (Triangle_ID_Only, Triangle_ID_Only);
};


enum DISABLING_POINTS_METHOD {
    NO_DISABLED_POINTS,
    DISABLE_POINTS_BY_INDEX,
    DISABLE_POINTS_BY_RANGE,
};


class Grid_info_manager {
private:
    double *coord_values[2];
    int num_points;
    double min_lon;
    double max_lon;
    double min_lat;
    double max_lat;
    bool is_cyclic;
    DISABLING_POINTS_METHOD disabling_method;
    int disabling_num;
    void* disabling_data;

    void gen_basic_grid();

public:
    /* for unittest */
    Grid_info_manager();
    virtual ~Grid_info_manager();
    virtual double** get_grid_coord_values(int);
    virtual bool* get_grid_mask(int);
    virtual int get_grid_num_points(int);
    virtual void get_grid_boundry(int, double*, double*, double*, double*);
    virtual void set_grid_boundry(int, double, double, double, double);
    virtual bool is_grid_cyclic(int);
    virtual bool read_grid_from_text(const char []);
    virtual void get_disabled_points_info(int, DISABLING_POINTS_METHOD*, int*, void**);
#ifdef NETCDF
    virtual void read_grid_from_nc(const char [], const char [], const char []);
    void gen_three_polar_grid();
    void gen_latlon_grid();
    void gen_latlon_90_grid();
#endif
};

extern Grid_info_manager *grid_info_mgr;
#endif
