/***************************************************************
  *  Copyright (c) 2013, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Haoyu Yang,
  *  supervised by Dr. Li Liu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "delaunay_grid_decomposition_mgt.h"
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

#define DEFAULT_NUM_POINTS_THRESHOLD 100
#define DEFAULT_EXPANGDING_RATIO 0.1
#define PDLN_TOLERABLE_ERROR 0.0001
#define PDLN_MAX_ITER_COUNT 10
#define PDLN_FLOAT_EQ_ERROR 1e-10

#define PDLN_DECOMPOSE_COMMON_MODE 0
#define PDLN_DECOMPOSE_SPOLAR_MODE 1
#define PDLN_DECOMPOSE_NPOLAR_MODE 2

#define PDLN_SPOLAR_MAX_LAT -45.0
#define PDLN_NPOLAR_MIN_LAT 45.0


Boundry& Boundry::operator= (Boundry& boundry) {
    min_lat = boundry.min_lat;
    min_lon = boundry.min_lon;
    max_lat = boundry.max_lat;
    max_lon = boundry.max_lon;
}


Search_tree_node::Search_tree_node(Search_tree_node *parent, double *coord_value[2], int num_points, Boundry boundry) {
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

    this->local_cells_global_index = new int[num_points];
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


/* assumption: ids are already sorted. */
void Search_tree_node::update_processing_units_id(int* ids, int num) 
{
    this->processing_units_id.clear();
    for(int i = 0; i < num; i ++)
        this->processing_units_id.push_back(ids[i]);
}


void Search_tree_node::update_processing_units_id(vector<int> proc_units_id) 
{
    this->processing_units_id.clear();
    this->processing_units_id = proc_units_id;
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


int Search_tree_node::do_decompose(Workload_info *workload_info, double *child_cells_coord[4], int child_num_cells[2],
                                   Boundry child_boundry[2], vector<int> child_proc_units_id[2], int mode)
{
    double length[2], boundry_values[4], child_total_workload[2];
    Midline midline;
    int i, child_cells_global_index[2];
    int iteration_count;

    assert(this->processing_units_id.size() > 1);
    boundry_values[0] = this->kernel_boundry->min_lon;
    boundry_values[1] = this->kernel_boundry->min_lat;
    boundry_values[2] = this->kernel_boundry->max_lon;
    boundry_values[3] = this->kernel_boundry->max_lat;
    assert(boundry_values[0] < boundry_values[2]);
    assert(boundry_values[1] < boundry_values[3]);
    length[0] = boundry_values[2] - boundry_values[0];
    length[1] = boundry_values[3] - boundry_values[1];
    if(length[0] < 0.0)
        length[0] += 360.0;
    assert(length[0] <= 360 && length[0] >= 0 && length[1] <= 180 && length[1] >= 0);
    
    if(mode == PDLN_DECOMPOSE_SPOLAR_MODE || mode == PDLN_DECOMPOSE_NPOLAR_MODE)
        midline.type = PDLN_LAT;
    else if(length[1] > length[0])
        midline.type = PDLN_LAT;
    else
        midline.type = PDLN_LON;

    child_proc_units_id[0].clear();
    child_proc_units_id[1].clear();
    if(mode == PDLN_DECOMPOSE_COMMON_MODE) {
        for(i = 0; i < this->processing_units_id.size()/2; i++)
            child_proc_units_id[0].push_back(this->processing_units_id[i]);
        for(; i < this->processing_units_id.size(); i++)
            child_proc_units_id[1].push_back(this->processing_units_id[i]);
    }
    else if(mode == PDLN_DECOMPOSE_SPOLAR_MODE) {
        child_proc_units_id[0].push_back(this->processing_units_id[0]);
        for(i = 1; i < this->processing_units_id.size(); i++)
            child_proc_units_id[1].push_back(this->processing_units_id[i]);
    }
    else if(mode == PDLN_DECOMPOSE_NPOLAR_MODE) {
        for(i = 0; i < this->processing_units_id.size()-1; i++)
            child_proc_units_id[0].push_back(this->processing_units_id[i]);
        child_proc_units_id[1].push_back(this->processing_units_id[i]);
    }
    else
        exit(2);

    assert(child_proc_units_id[0].size() + child_proc_units_id[1].size() == this->processing_units_id.size());
    if(this->processing_units_id.size() > 1) {
        for(i = 0, child_total_workload[0] = 0.0; i < child_proc_units_id[0].size(); i++)
            child_total_workload[0] += workload_info->workloads[child_proc_units_id[0][i]];
        for(i = 0, child_total_workload[1] = 0.0; i < child_proc_units_id[1].size(); i++)
            child_total_workload[1] += workload_info->workloads[child_proc_units_id[1][i]];
        //child_total_workload[0] = this->num_local_kernel_cells * child_proc_units_id[0].size() / this->processing_units_id.size();
        //child_total_workload[1] = this->num_local_kernel_cells * child_proc_units_id[1].size() / this->processing_units_id.size();

        //midline.value = boundry_values[midline.type] + length[midline.type] * child_proc_units_id[0].size() / this->processing_units_id.size();
        midline.value = boundry_values[midline.type] + length[midline.type] * child_total_workload[0] / (child_total_workload[0] + child_total_workload[1]);
        this->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);

        iteration_count = 1;
        if(child_num_cells[0] == 0)
            child_num_cells[0] = 1;
        if(child_num_cells[1] == 0)
            child_num_cells[1] = 1;
        while(fabs(child_num_cells[0]/child_num_cells[1] - child_total_workload[0]/child_total_workload[1]) > PDLN_TOLERABLE_ERROR) {
        //while(fabs(child_num_cells[0]/child_num_cells[1] - child_proc_units_id[0].size()/child_proc_units_id[1].size()) > PDLN_TOLERABLE_ERROR) {
            if(iteration_count++ > PDLN_MAX_ITER_COUNT)
                break;

            if(child_num_cells[0] < child_total_workload[0]) {
#ifdef DEBUG
                assert(child_num_cells[1] >= child_total_workload[1] || fabs(child_total_workload[1] - child_num_cells[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
                midline.value += (boundry_values[2+midline.type] - midline.value) * (child_num_cells[1] - child_total_workload[1]) / child_num_cells[1];
            }
            else {
#ifdef DEBUG
                assert(child_num_cells[1] <= child_total_workload[1] || fabs(child_total_workload[1] - child_num_cells[1]) < PDLN_FLOAT_EQ_ERROR);
#endif
                midline.value -= (midline.value - boundry_values[midline.type]) * (child_num_cells[0] - child_total_workload[0]) / child_num_cells[0];
            }
            assert(midline.value > boundry_values[midline.type]);
            assert(midline.value < boundry_values[2+midline.type]);
            /* TODO: Search only half of the whole points, but not the whole points */
            this->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);
            if(child_num_cells[0] == 0)
                child_num_cells[0] = 1;
            if(child_num_cells[1] == 0)
                child_num_cells[1] = 1;
        }
    }
    else
        midline.value = boundry_values[2+midline.type];
    
    if(midline.type == PDLN_LON) {
        child_boundry[0].min_lat = child_boundry[1].min_lat = this->kernel_boundry->min_lat;
        child_boundry[0].max_lat = child_boundry[1].max_lat = this->kernel_boundry->max_lat;
        child_boundry[0].min_lon = this->kernel_boundry->min_lon;
        child_boundry[0].max_lon = child_boundry[1].min_lon = midline.value;
        child_boundry[1].max_lon = this->kernel_boundry->max_lon;
    }
    else if(midline.type == PDLN_LAT) {
        child_boundry[0].min_lon = child_boundry[1].min_lon = this->kernel_boundry->min_lon;
        child_boundry[0].max_lon = child_boundry[1].max_lon = this->kernel_boundry->max_lon;
        child_boundry[0].min_lat = this->kernel_boundry->min_lat;
        child_boundry[0].max_lat = child_boundry[1].min_lat = midline.value;
        child_boundry[1].max_lat = this->kernel_boundry->max_lat;
    }
    else
        assert(false);

}


/* proc_info = processing_info_mgr->get_processing_info(component_id) */
Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_resource *proc_info)
{
    double **coord_values;
    Boundry boundry;
    int num_points;

    this->original_grid = grid_id;
    this->processing_info = proc_info;
    assert(this->processing_info != NULL);
    this->workload_info = NULL;

    this->min_num_points_per_chunk = DEFAULT_NUM_POINTS_THRESHOLD;
    
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);
    grid_info_mgr->get_grid_boundry(grid_id, &boundry.min_lat, &boundry.max_lat, &boundry.min_lon, &boundry.max_lon);
    this->processing_info->get_num_total_processing_units();
    this->search_tree_root = new Search_tree_node(NULL, coord_values, num_points, boundry);
    //this->initialze_workload();
}


int Delaunay_grid_decomposition::initialze_workload()
{
    int max_num_processing_units, num_actived_processing_units;
    int *actived_units_id;
    double average_workload;

    assert(this->min_num_points_per_chunk > 0);
    max_num_processing_units = (grid_info_mgr->get_grid_num_points(this->original_grid) + this->min_num_points_per_chunk - 1) / this->min_num_points_per_chunk;

    num_actived_processing_units = min(this->processing_info->get_num_total_processing_units(), max_num_processing_units);
    average_workload = (double)grid_info_mgr->get_grid_num_points(this->original_grid) / num_actived_processing_units;

    this->processing_info->pick_out_actived_processing_units(this->original_grid, num_actived_processing_units, average_workload);
    this->workload_info = this->processing_info->search_grid_workload_info(this->original_grid);
    this->search_tree_root->update_processing_units_id(this->workload_info->actived_common_id,
                                                       this->workload_info->size_actived);
}


/* "common_node" means non-polar node */
int Delaunay_grid_decomposition::decompose_common_node_recursively(Search_tree_node *node)
{
    double *child_cells_coord[4];
    int child_num_cells[2];
    Boundry child_boundry[2];
    vector<int> child_proc_id[2];

    assert(node->processing_units_id.size() > 0);
    if(node->processing_units_id.size() == 1) {
        if(this->have_local_processing_units_id(node->processing_units_id)) {
            this->local_leaf_nodes.push_back(node);
        }
        return 0;
    }
    
    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[node->num_local_kernel_cells];
    
    node->do_decompose(this->workload_info, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);

    node->first_child = new Search_tree_node(node, child_cells_coord,   child_num_cells[0], child_boundry[0]);
    node->third_child = new Search_tree_node(node, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

    /* child_proc_id[0] can be modified by this->workload_info->update_workloads */
    this->workload_info->update_workloads(child_num_cells[0], child_proc_id[0], this->min_num_points_per_chunk);
    this->workload_info->update_workloads(child_num_cells[1], child_proc_id[1], this->min_num_points_per_chunk);
    node->first_child->update_processing_units_id(child_proc_id[0]);
    node->third_child->update_processing_units_id(child_proc_id[1]);

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
    //TODO: optimize new delete
    
    if(this->have_local_processing_units_id(node->first_child->processing_units_id)) 
        this->decompose_common_node_recursively(node->first_child);
    if(this->have_local_processing_units_id(node->third_child->processing_units_id)) 
        this->decompose_common_node_recursively(node->third_child);

}


int Delaunay_grid_decomposition::assign_polars(bool assign_south_polar, bool assign_north_polar)
{
    double *child_cells_coord[4];
    Boundry child_boundry[2];
    int child_num_cells[2];
    vector<int> child_proc_id[2];
    Midline midline;
    
    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[this->search_tree_root->num_local_kernel_cells];

    if(assign_south_polar) {
        this->search_tree_root->do_decompose(this->workload_info, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_SPOLAR_MODE);
        if(child_boundry[0].max_lat > PDLN_SPOLAR_MAX_LAT || this->search_tree_root->processing_units_id.size() == 1) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_SPOLAR_MAX_LAT;
            this->search_tree_root->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);

            if(child_num_cells[0] < this->min_num_points_per_chunk) {
                for(int i = 0; i < 4; i++)
                    delete[] child_cells_coord[i];
                return 1;
            }
            child_boundry[0].max_lat = PDLN_SPOLAR_MAX_LAT;
            child_boundry[1].min_lat = PDLN_SPOLAR_MAX_LAT;
            /* actually this->search_tree_root->processing_units_id[0] means 0 here. */
            this->workload_info->workloads[this->search_tree_root->processing_units_id[0]] -= child_num_cells[0];
            child_proc_id[1].insert(child_proc_id[1].begin(), this->search_tree_root->processing_units_id[0]);
        }
        this->search_tree_root->first_child  = new Search_tree_node(this->search_tree_root, child_cells_coord,   child_num_cells[0], child_boundry[0]);
        this->search_tree_root->second_child = new Search_tree_node(this->search_tree_root, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

        this->workload_info->update_workloads(child_num_cells[0], child_proc_id[0], this->min_num_points_per_chunk);
        this->workload_info->update_workloads(child_num_cells[1], child_proc_id[1], this->min_num_points_per_chunk);
        this->search_tree_root->first_child->update_processing_units_id(child_proc_id[0]);
        this->search_tree_root->second_child->update_processing_units_id(child_proc_id[1]);

        this->current_tree_node = this->search_tree_root->second_child;
        
        if(this->have_local_processing_units_id(this->search_tree_root->first_child->processing_units_id))
            this->local_leaf_nodes.push_back(this->search_tree_root->first_child);
    }
    
    if(assign_north_polar) {
        this->current_tree_node->do_decompose(this->workload_info, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_NPOLAR_MODE);
        if(child_boundry[1].min_lat < PDLN_NPOLAR_MIN_LAT || this->current_tree_node->processing_units_id.size() == 1) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_NPOLAR_MIN_LAT;
            this->current_tree_node->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);

            if(child_num_cells[1] < this->min_num_points_per_chunk) {
                for(int i = 0; i < 4; i++)
                    delete[] child_cells_coord[i];
                return 1;
            }
            child_boundry[0].max_lat = PDLN_NPOLAR_MIN_LAT;
            child_boundry[1].min_lat = PDLN_NPOLAR_MIN_LAT;
            this->workload_info->workloads[this->search_tree_root->processing_units_id.back()] -= child_num_cells[1];
            child_proc_id[0].push_back(this->search_tree_root->processing_units_id.back());
        }
        if(this->search_tree_root->second_child)
            delete this->search_tree_root->second_child;
        this->search_tree_root->second_child = new Search_tree_node(this->search_tree_root, child_cells_coord,   child_num_cells[0], child_boundry[0]);
        this->search_tree_root->third_child  = new Search_tree_node(this->search_tree_root, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

        this->workload_info->update_workloads(child_num_cells[0], child_proc_id[0], this->min_num_points_per_chunk);
        this->workload_info->update_workloads(child_num_cells[1], child_proc_id[1], this->min_num_points_per_chunk);
        this->search_tree_root->second_child->update_processing_units_id(child_proc_id[0]);
        this->search_tree_root->third_child->update_processing_units_id(child_proc_id[1]);

        this->current_tree_node = this->search_tree_root->second_child;

        if(this->have_local_processing_units_id(this->search_tree_root->third_child->processing_units_id))
            this->local_leaf_nodes.push_back(this->search_tree_root->third_child);
    }

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
}


int Delaunay_grid_decomposition::assign_cyclic_grid_for_single_unit()
{
    double *child_cells_coord[4];
    Boundry child_boundry[2];
    int child_num_cells[2];
    vector<int> child_proc_id[2];
    Midline midline;

    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[this->current_tree_node->num_local_kernel_cells];

    assert(this->current_tree_node->processing_units_id.size() == 1);
    midline.type = PDLN_LON;
    midline.value = 180.0;
    child_boundry[0] = child_boundry[1] = *this->current_tree_node->kernel_boundry;
    child_boundry[0].max_lon = child_boundry[1].min_lon = 180.0;
    child_proc_id[0] = child_proc_id[1] = this->current_tree_node->processing_units_id;
    this->current_tree_node->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);
    this->workload_info->workloads[this->current_tree_node->processing_units_id[0]] -= child_num_cells[0] + child_num_cells[1];

    this->current_tree_node->first_child = new Search_tree_node(this->current_tree_node, child_cells_coord,   child_num_cells[0], child_boundry[0]);
    this->current_tree_node->third_child = new Search_tree_node(this->current_tree_node, child_cells_coord+2, child_num_cells[1], child_boundry[1]);
    this->current_tree_node->first_child->update_processing_units_id(child_proc_id[0]);
    this->current_tree_node->third_child->update_processing_units_id(child_proc_id[1]);

    assert(this->have_local_processing_units_id(child_proc_id[0]));
    this->local_leaf_nodes.push_back(this->current_tree_node->first_child);
    this->local_leaf_nodes.push_back(this->current_tree_node->third_child);

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
}

//TODO: get faster
bool Delaunay_grid_decomposition::have_local_processing_units_id(vector<int> units_id)
{
    for(int i = 0; i < units_id.size(); i++)
        for(int j = 0; j < this->processing_info->get_num_local_proc_processing_units(); j++)
            if(units_id[i] == this->processing_info->get_local_proc_common_id()[j])
                return true;

    return false;
}

int Delaunay_grid_decomposition::generate_grid_decomposition()
{
    int num_south_polar, num_north_polar;
    double *child_cells_coord[4];
    Boundry child_boundry[2];
    int child_num_cells[2];
    vector<int> child_proc_id[2];
    int i;

    this->initialze_workload();
    this->current_tree_node = this->search_tree_root;
    num_south_polar = grid_info_mgr->get_polar_points(this->original_grid, 'S');
    num_north_polar = grid_info_mgr->get_polar_points(this->original_grid, 'N');
    this->assign_polars(num_south_polar > 2, num_north_polar > 2);
    //for(int i = 0; i < this->workload_info->size_actived; i++)
    //    printf("common_id: %d, workload: %lf\n", i, this->workload_info->workloads[i]);

    if(this->current_tree_node->processing_units_id.size() == 1 && grid_info_mgr->is_grid_cyclic(this->original_grid)) {
        this->assign_cyclic_grid_for_single_unit();
        return 0;
    }

    this->decompose_common_node_recursively(this->current_tree_node);
}


int Delaunay_grid_decomposition::rotate_grid()
{
    
}

double** Grid_info_manager::get_grid_coord_values(int grid_id)
{
}
int Grid_info_manager::get_grid_num_points(int grid_id)
{
}
void Grid_info_manager::get_grid_boundry(int grid_id, double* min_lat, double* max_lat, double* min_lon, double* max_lon)
{
}
int Grid_info_manager::get_polar_points(int grid_id, char polar)
{
}
bool Grid_info_manager::is_grid_cyclic(int grid_id)
{
}
