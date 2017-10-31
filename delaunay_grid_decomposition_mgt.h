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
#include "delaunay_triangulation_2D.h"
#include "processing_unit_mgt.h"

struct Midline {
    int type;
    double value;
};

struct Boundry {
    double min_lat;
    double max_lat;
    double min_lon;
    double max_lon;

    Boundry& operator= (Boundry& boundry);
};

class Search_tree_node {
private:
    friend class Delaunay_grid_decomposition;
    Search_tree_node *parent;
    Search_tree_node *first_child, *second_child, *third_child;
    Boundry *kernel_boundry;
    Boundry *expanded_boundry;
    Boundry *rotated_kernel_boundry;
    Boundry *rotated_expanded_boundry;
    double expanding_ratio;
    Midline midline;
    double *local_cells_coord[2]; //(clean after having children)
    int    *local_cells_global_index; //(kernel + expanded)
    int    num_local_kernel_cells;
    int    num_local_expanded_cells;
    
    vector <int> processing_units_id;
    vector <Search_tree_node *> neighbors;
    //XY_triangulation_class *triangulation;
public:    
    void decompose_iteratively(double*, double **, int*, Boundry*, vector<int>*, int);
    void decompose_with_certain_line(Midline, double**, int*);

    Search_tree_node(Search_tree_node*, double**, int, Boundry);
    ~Search_tree_node();
    void update_processing_units_id(int);
    void update_processing_units_id(vector<int>);

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
    //int *num_local_nodes_per_thread;
    int num_total_nodes;
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
public:
    Delaunay_grid_decomposition(int, Processing_resource*);
    Delaunay_grid_decomposition(int, Processing_resource*, int);
    ~Delaunay_grid_decomposition();
    int generate_grid_decomposition();
    int expand_boundry();
    int generate_delaunay_triangulizition();
    vector<Search_tree_node*> get_local_leaf_nodes() {return this->local_leaf_nodes; };
    int generate_trianglulation_for_local_decomp();
    int generate_trianglulation_for_whole_grid();
};


class Grid_info_manager {
private:
    double *coord_values[2];
    int num_points;
public:
    /* for unittest */
    Grid_info_manager();
    virtual ~Grid_info_manager();
    virtual double** get_grid_coord_values(int);
    virtual int get_grid_num_points(int);
    virtual void get_grid_boundry(int, double*, double*, double*, double*);
    virtual int get_polar_points(int, char);
    virtual bool is_grid_cyclic(int);
};

extern Grid_info_manager *grid_info_mgr;
#endif
