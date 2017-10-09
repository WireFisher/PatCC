/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#ifndef DELAUNAY_GRID_DECOMPOSITION_MGT_H
#define DELAUNAY_GRID_DECOMPOSITION_MGT_H

#define PDELAUNAY_LON 0
#define PDELAUNAY_LAT 1
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

struct Search_tree_node {
    struct Search_tree_node *parent;
    struct Search_tree_node *first_child, *second_child, *third_child;
    struct Boundry *kernel_boundry;
    struct Boundry *expanded_boundry;
    struct Boundry *rotated_kernel_boundry;
    struct Boundry *rotated_expanded_boundry;
    double expanding_ratio;
    Midline midline;
    double *local_cells_coord[2]; //(clean after having children)
    int    *local_cells_global_index; //(kernel + expanded)
    int    num_local_kernel_cells;
    int    num_local_expanded_cells;
    
    vector <int> processing_units_id;
    vector <Search_tree_node *> neighbors;
    //XY_triangulation_class *triangulation;

    Search_tree_node(Search_tree_node*, double**, Boundry, int);
    ~Search_tree_node();
    void update_processing_units_id(int*, int);
    int decompose_with_certain_line(Midline, double**, int*);
    int decompose_automatically(Workload_info*);
};

class Delaunay_grid_decomposition {
private:
    Processing_info *processing_info;
    //Remap_grid_class *original_grid;
    int original_grid;
    Search_tree_node *search_tree_root;
    Search_tree_node *current_tree_node;
    Search_tree_node *local_threads_node;
    int *num_local_nodes_per_thread;
    int num_total_nodes;
    int *local_units_id;
    int min_num_points_per_chunk;
    //double min_length_single_chunk;
    
    int rotate_grid();
    int initialze_workload();
    int assign_polars(bool, bool);
public:
    Delaunay_grid_decomposition(int, Processing_info*);
    ~Delaunay_grid_decomposition();
    int generate_grid_decomposition();
    int expand_boundry();
    int generate_delaunay_triangulizition();
};


#endif
