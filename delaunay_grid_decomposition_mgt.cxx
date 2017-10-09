/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "delaunay_grid_decomposition_mgt.h"
#include "testcase.h"
#include <cstddef>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

#define DEFAULT_NUM_POINTS_THRESHOLD 100
#define DEFAULT_EXPANGDING_RATIO 0.1
#define PDELAUNAY_TOLERABLE_ERROR 0.0001
#define PDELAUNAY_MAX_ITER_COUNT 10

Boundry& Boundry::operator= (Boundry& boundry) {
    min_lat = boundry.min_lat;
    min_lon = boundry.min_lon;
    max_lat = boundry.max_lat;
    max_lon = boundry.max_lon;
}


Search_tree_node::Search_tree_node(Search_tree_node *parent, double *coord_value[2], Boundry boundry, int num_points) {
    this->parent = parent;
    this->first_child = this->second_child = this->third_child = NULL;

    this->kernel_boundry = new Boundry();
    this->expanded_boundry = new Boundry();
    *this->kernel_boundry = *this->expanded_boundry = boundry;
    this->rotated_kernel_boundry = this->rotated_expanded_boundry = NULL;

    this->expanding_ratio = DEFAULT_EXPANGDING_RATIO;
    this->midline.type = -1;
    this->midline.value = -361.0;
    this->local_cells_coord[0] = new double[num_points];
    this->local_cells_coord[1] = new double[num_points];
    memcpy(this->local_cells_coord[0], coord_value[0], num_points * sizeof(double));
    memcpy(this->local_cells_coord[1], coord_value[1], num_points * sizeof(double));

    for(int i=0; i < num_points; i++)
        this->local_cells_global_index[i] = i;

    this->num_local_kernel_cells = this->num_local_expanded_cells = num_points;
    //this->triangulation = NULL;
}


Search_tree_node::~Search_tree_node()
{
    delete[] this->local_cells_coord[0];
    delete[] this->local_cells_coord[1];
    delete this->kernel_boundry;
    delete this->expanded_boundry;
    if(this->rotated_kernel_boundry)
        delete this->rotated_kernel_boundry;
    if(this->rotated_expanded_boundry)
        delete this->rotated_expanded_boundry;
}


void Search_tree_node::update_processing_units_id(int* ids,int num) 
{
    this->processing_units_id.clear();
    for(int i; i < num; i ++)
        this->processing_units_id.push_back(ids[i]);
}

/* proc_info = processing_info_mgr->get_processing_info(component_id) */
Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_info *proc_info)
{
    double **coord_values;
    Boundry boundry;
    int num_points;

    this->original_grid = grid_id;
    this->processing_info = proc_info;
    //assert this->processing_info != NULL

    this->min_num_points_per_chunk = DEFAULT_NUM_POINTS_THRESHOLD;
    
    coord_values = get_grid_coord_values(grid_id);
    num_points = get_grid_num_points(grid_id);
    get_grid_boundry(grid_id, &boundry.min_lat, &boundry.max_lat, &boundry.min_lon, &boundry.max_lon);
    this->processing_info->get_num_total_processing_units();
    this->search_tree_root = new Search_tree_node(NULL, coord_values, boundry, num_points);
    //this->initialze_workload();
}


int Delaunay_grid_decomposition::initialze_workload()
{
    int max_num_processing_units, num_actived_processing_units;
    int *actived_units_id;
    double average_workload;

    //assert min_num_points_per_chunk && min_num_points_per_chunk > 0
    max_num_processing_units = (get_grid_size(this->original_grid) + this->min_num_points_per_chunk - 1) / this->min_num_points_per_chunk;

    num_actived_processing_units = min(this->processing_info->get_num_total_processing_units(), max_num_processing_units);
    average_workload = (double)get_grid_size(this->original_grid) / num_actived_processing_units;

    this->processing_info->pick_out_actived_processing_units(this->original_grid, num_actived_processing_units, average_workload);
    this->search_tree_root->update_processing_units_id(this->processing_info->search_grid_workload_info(this->original_grid)->actived_common_id,
                                                       this->processing_info->search_grid_workload_info(this->original_grid)->size_actived);
}


int Search_tree_node::decompose_with_certain_line(Midline midline, double *child_cells_coord[4], int child_num_cells[2])
{
    child_num_cells[0] = 0;
    child_num_cells[1] = 0;
    for(int i = 0; i < this->num_local_kernel_cells; i ++) {
        if(local_cells_coord[midline.type][i] < midline.value) {
            child_cells_coord[midline.type][child_num_cells[0]] = local_cells_coord[midline.type][i];
            child_cells_coord[(midline.type+1)%2][child_num_cells[0]++] = local_cells_coord[(midline.type+1)%2][i];
        }
        else {
            child_cells_coord[2+midline.type][child_num_cells[1]] = local_cells_coord[midline.type][i];
            child_cells_coord[2+(midline.type+1)%2][child_num_cells[1]++] = local_cells_coord[(midline.type+1)%2][i];
        }
    }
}


int Search_tree_node::decompose_automatically(Workload_info *workload_info)
{
    double length[2], boundry_values[4], child_total_workload[2], *child_cells_coord[4];
    Midline midline;
    vector<int> child_proc_units_id[2];
    int i, child_cells_global_index[2], child_num_cells[2];
    int iteration_count;
    Boundry child_boundry[2];

    if(this->processing_units_id.size() == 1)
        return

    //assert this->processing_units_id.size() > 1
    boundry_values[0] = this->kernel_boundry->min_lon;
    boundry_values[1] = this->kernel_boundry->min_lat;
    boundry_values[2] = this->kernel_boundry->max_lon;
    boundry_values[3] = this->kernel_boundry->max_lat;
    length[0] = boundry_values[2] - boundry_values[0];
    length[1] = boundry_values[3] - boundry_values[1];
    if(length[0] < 0.0)
        length[0] += 360.0;
    //assert length_lon <= 360 && length_lon >= 0 && length_lat <= 90 && length_lat >= -90
    if(length[1] > length[0])
        midline.type = PDELAUNAY_LAT;
    else
        midline.type = PDELAUNAY_LON;

    for(i = 0; i < this->processing_units_id.size()/2; i++)
        child_proc_units_id[0].push_back(this->processing_units_id[i]);
    for(; i < this->processing_units_id.size(); i++)
        child_proc_units_id[1].push_back(this->processing_units_id[i]);
    //assert child_proc_units_id[0].size + child_proc_units_id[1].size == this->processing_units_id.size()
    //assert child_proc_units_id[0].size <= child_proc_units_id[1].size
    midline.value = boundry_values[midline.type] + length[midline.type] * child_proc_units_id[0].size() / this->processing_units_id.size();

    for(i = 0; i < 4; i++)
        child_cells_coord[i] = new double[this->num_local_kernel_cells];
    for(i = 0, child_total_workload[0] = 0.0; i < child_proc_units_id[0].size(); i++)
        child_total_workload[0] += workload_info->workloads[child_proc_units_id[0][i]];
    for(i = 0, child_total_workload[1] = 0.0; i < child_proc_units_id[1].size(); i++)
        child_total_workload[1] += workload_info->workloads[child_proc_units_id[1][i]];
    //child_total_workload[0] = this->num_local_kernel_cells * child_proc_units_id[0].size() / this->processing_units_id.size();
    //child_total_workload[1] = this->num_local_kernel_cells * child_proc_units_id[1].size() / this->processing_units_id.size();

    this->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);

    iteration_count = 1;
    while(fabs(child_num_cells[0]/child_num_cells[1] - child_proc_units_id[0].size()/child_proc_units_id[1].size()) > PDELAUNAY_TOLERABLE_ERROR) {
        if(iteration_count++ > PDELAUNAY_MAX_ITER_COUNT)
            break;

        if(child_num_cells[0] < child_total_workload[0]) {
            //assert child_num_cells[1] > child_total_workload[1]
            midline.value += (boundry_values[2+midline.type] - midline.value) * (child_num_cells[1] - child_total_workload[1]) / child_num_cells[1];
        }
        else {
            //assert child_num_cells[1] <= child_total_workload[1]
            midline.value -= (midline.value - boundry_values[midline.type]) * (child_num_cells[0] - child_total_workload[0]) / child_num_cells[0];
        }
        /* TODO: Search only half of the whole points, but not the whole points */
        this->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);
    }
    
    if(midline.type == PDELAUNAY_LON) {
        child_boundry[0].min_lat = child_boundry[1].min_lat = this->kernel_boundry->min_lat;
        child_boundry[0].max_lat = child_boundry[1].max_lat = this->kernel_boundry->max_lat;
        child_boundry[0].min_lon = this->kernel_boundry->min_lon;
        child_boundry[0].max_lon = child_boundry[1].min_lon = midline.value;
        child_boundry[1].max_lon = this->kernel_boundry->max_lon;
    }
    else if(midline.type == PDELAUNAY_LAT) {
        child_boundry[0].min_lon = child_boundry[1].min_lon = this->kernel_boundry->min_lon;
        child_boundry[0].max_lon = child_boundry[1].max_lon = this->kernel_boundry->max_lon;
        child_boundry[0].min_lat = this->kernel_boundry->min_lat;
        child_boundry[0].max_lat = child_boundry[1].min_lat = midline.value;
        child_boundry[1].max_lat = this->kernel_boundry->max_lat;
    }
    else
        //assert false
        exit(1);
    this->first_child = new Search_tree_node(this, child_cells_coord,   child_boundry[0], child_num_cells[0]);
    this->third_child = new Search_tree_node(this, child_cells_coord+2, child_boundry[1], child_num_cells[1]);
    //update units id

    for(i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
}


int Delaunay_grid_decomposition::assign_polars(bool assign_south_polar, bool assign_north_polar)
{
    
}


int Delaunay_grid_decomposition::generate_grid_decomposition()
{
    int num_south_polar, num_north_polar;

    this->initialze_workload();
    this->current_tree_node = this->search_tree_root;
    num_south_polar = get_polar_points('S');
    num_north_polar = get_polar_points('N');
    //this->assign_polars();
    
}


int Delaunay_grid_decomposition::rotate_grid()
{
    
}
