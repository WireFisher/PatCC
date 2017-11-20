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

#define DEFAULT_EXPANGDING_RATIO 0.1
#define PDLN_EXPECTED_EXPANDING_LOOP_TIMES 3

#define PDLN_SPOLAR_MAX_LAT -45.0
#define PDLN_NPOLAR_MIN_LAT 45.0

#define PDLN_MAX_ITER_COUNT 10

#define PDLN_TOLERABLE_ERROR 0.0001
#define PDLN_FLOAT_EQ_ERROR 1e-10

#define PDLN_DECOMPOSE_COMMON_MODE 0
#define PDLN_DECOMPOSE_SPOLAR_MODE 1
#define PDLN_DECOMPOSE_NPOLAR_MODE 2

#define PDLN_NODE_TYPE_COMMON PDLN_DECOMPOSE_COMMON_MODE
#define PDLN_NODE_TYPE_SPOLAR PDLN_DECOMPOSE_SPOLAR_MODE
#define PDLN_NODE_TYPE_NPOLAR PDLN_DECOMPOSE_NPOLAR_MODE


Boundry& Boundry::operator= (Boundry& boundry)
{
    min_lat = boundry.min_lat;
    min_lon = boundry.min_lon;
    max_lat = boundry.max_lat;
    max_lon = boundry.max_lon;
    return *this;
}


Boundry& Boundry::operator* (double ratio)
{
    min_lat *= ratio;
    min_lon *= ratio;
    max_lat *= ratio;
    max_lon *= ratio;
    return *this;
}


void Boundry::legalize()
{
    min_lat = std::max(min_lat, -90.0);
    min_lon = std::max(min_lon, 0.0);
    max_lat = std::min(max_lat, 90.0);
    max_lon = std::min(max_lon, 360.0);
}


void Boundry::legalize(Boundry outer_boundry)
{
    min_lat = std::max(min_lat, outer_boundry.min_lat);
    min_lon = std::max(min_lon, outer_boundry.min_lon);
    max_lat = std::min(max_lat, outer_boundry.max_lat);
    max_lon = std::min(max_lon, outer_boundry.max_lon);
}


Search_tree_node::Search_tree_node(Search_tree_node *parent, double *coord_value[2], int num_points, Boundry boundry) {
    this->parent = parent;
    this->first_child = this->second_child = this->third_child = NULL;

    this->kernel_boundry = new Boundry();
    this->expanded_boundry = new Boundry();
    *this->kernel_boundry = *this->expanded_boundry = boundry;
    this->rotated_kernel_boundry = this->rotated_expanded_boundry = NULL;

    //this->expanding_ratio = DEFAULT_EXPANGDING_RATIO;
    this->midline.type = -1;
    this->midline.value = -361.0;
    this->local_cells_coord[0] = new double[num_points];
    this->local_cells_coord[1] = new double[num_points];
    memcpy(this->local_cells_coord[0], coord_value[0], num_points * sizeof(double));
    memcpy(this->local_cells_coord[1], coord_value[1], num_points * sizeof(double));
    this->expanded_cells_coord[0] = this->expanded_cells_coord[1] = NULL;
    this->len_expanded_cells_coord_buf = 0;

    /*
    this->local_cells_global_index = new int[num_points];
    for(int i=0; i < num_points; i++)
        this->local_cells_global_index[i] = i;
    */

    this->num_local_kernel_cells = num_points;
    this->num_local_expanded_cells = 0;
    this->node_type = PDLN_NODE_TYPE_COMMON;
    this->triangulation = NULL;
}


Search_tree_node::~Search_tree_node()
{
    delete this->kernel_boundry;
    delete this->expanded_boundry;
    delete[] this->local_cells_coord[0];
    delete[] this->local_cells_coord[1];
    if(this->rotated_kernel_boundry != NULL)
        delete this->rotated_kernel_boundry;
    if(this->rotated_expanded_boundry != NULL)
        delete this->rotated_expanded_boundry;
    //delete[] this->local_cells_global_index;
    if(this->first_child != NULL)
        delete this->first_child;
    if(this->second_child != NULL)
        delete this->second_child;
    if(this->third_child != NULL)
        delete this->third_child;
}


/* assumption: ids are already sorted. */
/*
void Search_tree_node::update_processing_units_id(int* ids, int num) 
{
    this->processing_units_id.clear();
    for(int i = 0; i < num; i ++)
        this->processing_units_id.push_back(ids[i]);
}
*/

void Search_tree_node::update_processing_units_id(int num) 
{
    this->processing_units_id.clear();
    for(int i = 0; i < num; i ++)
        this->processing_units_id.push_back(i);
}


void Search_tree_node::update_processing_units_id(vector<int> proc_units_id) 
{
    this->processing_units_id.clear();
    this->processing_units_id = proc_units_id;
}


void Search_tree_node::generate_local_triangulation()
{
    std::vector<Point<double> > points;
    //for(int i = 0; i < num_local_kernel_cells; i++)
        //points.push_back(Vector2<double>(local_cells_coord[0][i], local_cells_coord[1][i]));
    //for(int i = 0; i < num_local_expanded_cells; i++)
        //points.push_back(Vector2<double>(expanded_cells_coord[0][i], expanded_cells_coord[1][i]));

    this->triangulation = new Delaunay<double>();
    this->triangulation->triangulate(points);
}

void Search_tree_node::decompose_with_certain_line(Midline midline, double *child_cells_coord[4], int child_num_cells[2])
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


void Search_tree_node::decompose_iteratively(double *workloads, double *child_cells_coord[4], int child_num_cells[2],
                                   Boundry child_boundry[2], vector<int> child_proc_units_id[2], int mode)
{
    double length[2], boundry_values[4], child_total_workload[2];
    Midline midline;
    unsigned int i;
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

    /* PDLN_DECOMPOSE_COMMON_MODE: 0   1   2   3 | 4   5   6   7
     * PDLN_DECOMPOSE_SPOLAR_MODE: 0 | 1   2   3   4   5   6   7
     * PDLN_DECOMPOSE_NPOLAR_MODE: 0   1   2   3   4   5   6 | 7 */
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
        assert(false);

    assert(child_proc_units_id[0].size() + child_proc_units_id[1].size() == this->processing_units_id.size());
    if(this->processing_units_id.size() > 1) {
        for(i = 0, child_total_workload[0] = 0.0; i < child_proc_units_id[0].size(); i++)
            child_total_workload[0] += workloads[child_proc_units_id[0][i]];
        for(i = 0, child_total_workload[1] = 0.0; i < child_proc_units_id[1].size(); i++)
            child_total_workload[1] += workloads[child_proc_units_id[1][i]];

        midline.value = boundry_values[midline.type] + length[midline.type] * child_total_workload[0] / (child_total_workload[0] + child_total_workload[1]);
        this->decompose_with_certain_line(midline, child_cells_coord, child_num_cells);
        assert(child_num_cells[0] != 0 && child_num_cells[1] != 0);

        iteration_count = 1;
        while(fabs(child_num_cells[0]/child_num_cells[1] - child_total_workload[0]/child_total_workload[1]) > PDLN_TOLERABLE_ERROR) {
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
            assert(child_num_cells[0] != 0 && child_num_cells[1] != 0);
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


void Search_tree_node::add_expanded_points(double *coord_value[2], int num_points)
{
    double *tmp_coord_value[2];
    if(num_local_expanded_cells + num_points > len_expanded_cells_coord_buf) {
        len_expanded_cells_coord_buf = num_local_expanded_cells + num_points * 4 * PDLN_EXPECTED_EXPANDING_LOOP_TIMES;
        tmp_coord_value[0] = new double[len_expanded_cells_coord_buf];
        tmp_coord_value[1] = new double[len_expanded_cells_coord_buf];
        memcpy(tmp_coord_value[0], local_cells_coord[0], sizeof(double) * num_local_expanded_cells);
        memcpy(tmp_coord_value[1], local_cells_coord[1], sizeof(double) * num_local_expanded_cells);
        delete[] local_cells_coord[0];
        delete[] local_cells_coord[1];
        local_cells_coord[0] = tmp_coord_value[0];
        local_cells_coord[1] = tmp_coord_value[1];
    }
    memcpy(local_cells_coord[0] + num_local_expanded_cells, coord_value[0], sizeof(double) * num_points);
    memcpy(local_cells_coord[1] + num_local_expanded_cells, coord_value[1], sizeof(double) * num_points);
    num_local_expanded_cells += num_points;
    assert(num_local_expanded_cells >= num_points);
}


void Search_tree_node::add_neighbors(vector<Search_tree_node*> neighbors)
{
    this->neighbors = neighbors;
    for(unsigned int i = 0; i < neighbors.size(); i++)
        if(find(this->neighbors.begin(), this->neighbors.end(), neighbors[i]) == this->neighbors.end())
            this->neighbors.push_back(neighbors[i]);
}


bool Search_tree_node::check_expanded_triangle_consistency()
{
    if(this->triangulation == NULL)
        return false;
}


Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_resource *proc_info, int min_num_points_per_chunk)
{
    double **coord_values;
    Boundry boundry;
    int num_points;

    this->original_grid = grid_id;
    this->processing_info = proc_info;
    assert(this->processing_info != NULL);

    this->min_num_points_per_chunk = min_num_points_per_chunk;
    
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);
    grid_info_mgr->get_grid_boundry(grid_id, &boundry.min_lat, &boundry.max_lat, &boundry.min_lon, &boundry.max_lon);
    this->processing_info->get_num_total_processing_units();
    this->search_tree_root = new Search_tree_node(NULL, coord_values, num_points, boundry);
    this->active_processing_common_id = NULL;
    this->workloads = NULL;
    //this->initialze_workload();
}

/* TODO: Delete this */
Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_resource *proc_info)
{
    double **coord_values;
    Boundry boundry;
    int num_points;

    this->original_grid = grid_id;
    this->processing_info = proc_info;
    assert(this->processing_info != NULL);

    this->min_num_points_per_chunk = 100;
    
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);
    grid_info_mgr->get_grid_boundry(grid_id, &boundry.min_lat, &boundry.max_lat, &boundry.min_lon, &boundry.max_lon);
    this->processing_info->get_num_total_processing_units();
    this->search_tree_root = new Search_tree_node(NULL, coord_values, num_points, boundry);
    this->active_processing_common_id = NULL;
    this->workloads = NULL;
    //this->initialze_workload();
}


Delaunay_grid_decomposition::~Delaunay_grid_decomposition()
{
    delete this->search_tree_root;
    delete[] this->active_processing_common_id;
    delete[] this->workloads; 
}


void Delaunay_grid_decomposition::initialze_workload()
{
    int max_num_processing_units, num_active_processing_units;
    int i, j;
    double average_workload;
    bool* active_units_flag;

    assert(this->min_num_points_per_chunk > 0);
    max_num_processing_units = (grid_info_mgr->get_grid_num_points(this->original_grid) + this->min_num_points_per_chunk - 1) / this->min_num_points_per_chunk;

    num_active_processing_units = std::min(this->processing_info->get_num_total_processing_units(), max_num_processing_units);
    average_workload = (double)grid_info_mgr->get_grid_num_points(this->original_grid) / num_active_processing_units;

    active_units_flag = new bool[this->processing_info->get_num_total_processing_units()];
    if(this->active_processing_common_id != NULL)
        delete[] this->active_processing_common_id;
    if(this->workloads != NULL)
        delete[] this->workloads;

    this->processing_info->pick_out_active_processing_units(num_active_processing_units, active_units_flag);

    this->active_processing_common_id = new int[num_active_processing_units];
    this->workloads = new double[num_active_processing_units];
    for(i = 0, j = 0; i < this->processing_info->get_num_total_processing_units(); i++)
        if(active_units_flag[i]) {
            this->workloads[j] = average_workload;
            this->active_processing_common_id[j++] = i;
        }

    assert(j == num_active_processing_units);

    this->search_tree_root->update_processing_units_id(num_active_processing_units);

    delete[] active_units_flag;
}


void Delaunay_grid_decomposition::update_workloads(int total_workload, vector<int> &ids)
{
    double old_total_workload = 0.0;
    double unassigned_workload = 0.0;

    for(unsigned int i = 0; i < ids.size(); i++)
        old_total_workload += this->workloads[ids[i]];

    for(unsigned int i = 0; i < ids.size(); i++)
        this->workloads[ids[i]] = this->workloads[ids[i]] * total_workload / old_total_workload;

    for(unsigned int i = 0; i < ids.size(); i++)
        if(this->workloads[ids[i]] < this->min_num_points_per_chunk) {
            unassigned_workload += this->workloads[ids[i]];
            ids.erase(ids.begin() + i);
        }

    assert(!0.0 > 0);
    if(unassigned_workload > 0)
        for(unsigned int i = 0; i < ids.size(); i++)
            this->workloads[ids[i]] += unassigned_workload / ids.size();
}


/* "common_node" means non-polar node */
void Delaunay_grid_decomposition::decompose_common_node_recursively(Search_tree_node *node)
{
    double *child_cells_coord[4];
    int child_num_cells[2];
    Boundry child_boundry[2];
    vector<int> child_proc_id[2];

    assert(node->processing_units_id.size() > 0);
    if(node->processing_units_id.size() == 1) {
        if(this->have_local_processing_units_id(node->processing_units_id))
            this->local_leaf_nodes.push_back(node);
        return;
    }
    
    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[node->num_local_kernel_cells];
    
    node->decompose_iteratively(this->workloads, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);

    node->first_child = new Search_tree_node(node, child_cells_coord,   child_num_cells[0], child_boundry[0]);
    node->third_child = new Search_tree_node(node, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

    /* child_proc_id[0] can be modified by this->update_workloads */
    this->update_workloads(child_num_cells[0], child_proc_id[0]);
    this->update_workloads(child_num_cells[1], child_proc_id[1]);
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


bool Delaunay_grid_decomposition::two_regions_overlap(Boundry region1, Boundry region2)
{
    if(region1.max_lat <= region2.min_lat || region1.min_lat >= region2.max_lat)
        return false;
    if(region1.max_lon <= region2.min_lon || region1.min_lon >= region2.max_lon)
        return false;

    return true;
}


void Delaunay_grid_decomposition::search_leaf_nodes_overlapping_with_region_recursively(Search_tree_node *node, Boundry region, vector<Search_tree_node*> &leaf_nodes_found)
{
    double *child_cells_coord[4];
    int child_num_cells[2];
    Boundry child_boundry[2];
    vector<int> child_proc_id[2];

    assert(node->processing_units_id.size() > 0);
    if(node->processing_units_id.size() == 1) {
        if(this->two_regions_overlap(region, *node->kernel_boundry))
            leaf_nodes_found.push_back(node);
        return;
    }

    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[node->num_local_kernel_cells];
    
    node->decompose_iteratively(this->workloads, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);

    node->first_child = new Search_tree_node(node, child_cells_coord,   child_num_cells[0], child_boundry[0]);
    node->third_child = new Search_tree_node(node, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

    /* child_proc_id[0] can be modified by this->update_workloads */
    this->update_workloads(child_num_cells[0], child_proc_id[0]);
    this->update_workloads(child_num_cells[1], child_proc_id[1]);
    node->first_child->update_processing_units_id(child_proc_id[0]);
    node->third_child->update_processing_units_id(child_proc_id[1]);

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
    //TODO: optimize new delete

    if(this->two_regions_overlap(region, *node->first_child->kernel_boundry))
        this->search_leaf_nodes_overlapping_with_region_recursively(node->first_child, region, leaf_nodes_found);
    if(node->second_child != NULL)
        if(this->two_regions_overlap(region, *node->second_child->kernel_boundry))
            this->search_leaf_nodes_overlapping_with_region_recursively(node->second_child, region, leaf_nodes_found);
    if(this->two_regions_overlap(region, *node->third_child->kernel_boundry))
        this->search_leaf_nodes_overlapping_with_region_recursively(node->third_child, region, leaf_nodes_found);
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
        this->search_tree_root->decompose_iteratively(this->workloads, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_SPOLAR_MODE);
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
            this->workloads[this->search_tree_root->processing_units_id[0]] -= child_num_cells[0];
            child_proc_id[1].insert(child_proc_id[1].begin(), this->search_tree_root->processing_units_id[0]);
        }
        this->search_tree_root->first_child  = new Search_tree_node(this->search_tree_root, child_cells_coord,   child_num_cells[0], child_boundry[0]);
        this->search_tree_root->second_child = new Search_tree_node(this->search_tree_root, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

        this->search_tree_root->first_child->node_type = PDLN_NODE_TYPE_SPOLAR;

        this->update_workloads(child_num_cells[0], child_proc_id[0]);
        this->update_workloads(child_num_cells[1], child_proc_id[1]);
        this->search_tree_root->first_child->update_processing_units_id(child_proc_id[0]);
        this->search_tree_root->second_child->update_processing_units_id(child_proc_id[1]);

        this->current_tree_node = this->search_tree_root->second_child;
        
        if(this->have_local_processing_units_id(this->search_tree_root->first_child->processing_units_id))
            this->local_leaf_nodes.push_back(this->search_tree_root->first_child);
    }
    
    if(assign_north_polar) {
        this->current_tree_node->decompose_iteratively(this->workloads, child_cells_coord, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_NPOLAR_MODE);
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
            this->workloads[this->search_tree_root->processing_units_id.back()] -= child_num_cells[1];
            child_proc_id[0].push_back(this->search_tree_root->processing_units_id.back());
        }
        if(this->search_tree_root->second_child != NULL)
            delete this->search_tree_root->second_child;
        this->search_tree_root->second_child = new Search_tree_node(this->search_tree_root, child_cells_coord,   child_num_cells[0], child_boundry[0]);
        this->search_tree_root->third_child  = new Search_tree_node(this->search_tree_root, child_cells_coord+2, child_num_cells[1], child_boundry[1]);

        this->search_tree_root->third_child->node_type = PDLN_NODE_TYPE_NPOLAR;

        this->update_workloads(child_num_cells[0], child_proc_id[0]);
        this->update_workloads(child_num_cells[1], child_proc_id[1]);
        this->search_tree_root->second_child->update_processing_units_id(child_proc_id[0]);
        this->search_tree_root->third_child->update_processing_units_id(child_proc_id[1]);

        this->current_tree_node = this->search_tree_root->second_child;

        if(this->have_local_processing_units_id(this->search_tree_root->third_child->processing_units_id))
            this->local_leaf_nodes.push_back(this->search_tree_root->third_child);
    }

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
    return 0;
}


void Delaunay_grid_decomposition::assign_cyclic_grid_for_single_unit()
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
    this->workloads[this->current_tree_node->processing_units_id[0]] -= child_num_cells[0] + child_num_cells[1];

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
bool Delaunay_grid_decomposition::have_local_processing_units_id(vector<int> chunk_id)
{
    for(unsigned int i = 0; i < chunk_id.size(); i++)
        for(int j = 0; j < this->processing_info->get_num_local_proc_processing_units(); j++)
            if(this->active_processing_common_id[chunk_id[i]] == this->processing_info->get_local_proc_common_id()[j])
                return true;

    return false;
}

int Delaunay_grid_decomposition::generate_grid_decomposition()
{
    int num_south_polar, num_north_polar;

    this->initialze_workload();
    this->current_tree_node = this->search_tree_root;
    num_south_polar = grid_info_mgr->get_polar_points(this->original_grid, 'S');
    num_north_polar = grid_info_mgr->get_polar_points(this->original_grid, 'N');
    if(this->assign_polars(num_south_polar > 2, num_north_polar > 2))
        return 1;
    //for(int i = 0; i < this->workload_info->size_active; i++)
    //    printf("common_id: %d, workload: %lf\n", i, this->workload_info->workloads[i]);

    if(this->current_tree_node->processing_units_id.size() == 1 && grid_info_mgr->is_grid_cyclic(this->original_grid)) {
        this->assign_cyclic_grid_for_single_unit();
        return 0;
    }

    this->decompose_common_node_recursively(this->current_tree_node);
    return 0;
}


int Delaunay_grid_decomposition::rotate_grid()
{
    return 0;
}


void Delaunay_grid_decomposition::transform_into_rectangle(Boundry inner_boundry, Boundry outer_boundry, Boundry sub_rectangle[4])
{
    sub_rectangle[0].min_lat = sub_rectangle[1].min_lat = outer_boundry.min_lat;
    sub_rectangle[1].max_lon = sub_rectangle[2].max_lon = outer_boundry.max_lon;
    sub_rectangle[2].max_lat = sub_rectangle[3].max_lat = outer_boundry.max_lat;
    sub_rectangle[3].min_lon = sub_rectangle[0].min_lon = outer_boundry.min_lon;

    sub_rectangle[0].max_lon = sub_rectangle[1].min_lon = inner_boundry.min_lon;
    sub_rectangle[1].max_lat = sub_rectangle[2].min_lat = inner_boundry.min_lat;
    sub_rectangle[2].min_lon = sub_rectangle[3].max_lon = inner_boundry.max_lon;
    sub_rectangle[3].min_lat = sub_rectangle[0].max_lat = inner_boundry.max_lat;
}


vector<Search_tree_node*> Delaunay_grid_decomposition::search_points_in_region(Boundry region, double *coord_values[2], int *num_points_found)
{
    vector<Search_tree_node*> leaf_nodes_found;
    
    *num_points_found = 0;
    if(fabs(region.min_lat - region.max_lat) < PDLN_FLOAT_EQ_ERROR || fabs(region.min_lon - region.max_lon) < PDLN_FLOAT_EQ_ERROR)
        return leaf_nodes_found;

    search_leaf_nodes_overlapping_with_region_recursively(this->search_tree_root, region, leaf_nodes_found);
    for(unsigned int i = 0; i < leaf_nodes_found.size(); i++)
        for(int j = 0; j < leaf_nodes_found[i]->num_local_kernel_cells; j++)
            if(leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] < region.max_lon &&
               leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] >= region.min_lon &&
               leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] < region.max_lat &&
               leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] >= region.min_lat) {
                coord_values[PDLN_LON][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j];
                coord_values[PDLN_LAT][(*num_points_found)++] = leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j];
            }
    return leaf_nodes_found;
}


int Delaunay_grid_decomposition::expand_tree_node_boundry(Search_tree_node* tree_node, double expanding_ratio)
{
    Boundry old_boundry = *tree_node->expanded_boundry;
    Boundry sub_rectangles[4];
    double *expanded_cells_coord[2];
    int num_points_found;

    *tree_node->expanded_boundry = *tree_node->expanded_boundry * expanding_ratio;
    //if(!is_syslic)
    tree_node->expanded_boundry->legalize(/* original grid boundry*/);

    if(tree_node->node_type == PDLN_NODE_TYPE_COMMON && old_boundry.max_lat >= 90 && tree_node->expanded_boundry->max_lat >= 90)
        return 1;
    if(tree_node->node_type == PDLN_NODE_TYPE_COMMON && old_boundry.min_lat <= -90 && tree_node->expanded_boundry->min_lat <= -90)
        return 1;

    transform_into_rectangle(old_boundry, *tree_node->expanded_boundry, sub_rectangles);

    expanded_cells_coord[0] = new double[search_tree_root->num_local_kernel_cells]; //FIXME: buf too large
    expanded_cells_coord[1] = new double[search_tree_root->num_local_kernel_cells];
    for(int i = 0; i < 4; i++){
        num_points_found = 0;
        vector<Search_tree_node*> leaf_nodes_found;
        leaf_nodes_found = search_points_in_region(sub_rectangles[i], expanded_cells_coord, &num_points_found);
        tree_node->add_expanded_points(expanded_cells_coord, num_points_found);
        tree_node->add_neighbors(leaf_nodes_found);
    }

    delete[] expanded_cells_coord[0];
    delete[] expanded_cells_coord[1];
    return 0;
}


int Delaunay_grid_decomposition::generate_trianglulation_for_local_decomp()
{
    //TODO: openmp parallel
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        int ret;
        double expanding_ratio = DEFAULT_EXPANGDING_RATIO;
        while(!local_leaf_nodes[i]->check_expanded_triangle_consistency()) {
            ret = expand_tree_node_boundry(local_leaf_nodes[i], expanding_ratio);
            //mpi bcast ret (NOTE: local threads, using tag)
            //if(one of the rets != 0) {
            //    expand_fail = true;
            //    break;
            //}
            if(local_leaf_nodes[i]->node_type) {
                //if(expanded_boundry not in correct region)
                //    return 2;
                //rotate
            }
            //do 2D delaunay triangulation
            local_leaf_nodes[i]->generate_local_triangulation();
            //expanding_ratio = DEFAULT_EXPANGDING_RATIO + 0.1
        }
    }
    //if(expand_fail)
    //    return 1;
    return 0;
}

int Delaunay_grid_decomposition::generate_trianglulation_for_whole_grid()
{
    return 0;
}

Grid_info_manager::Grid_info_manager()
{
    int size = 300;

    num_points = size * size;
    coord_values[0] = new double[num_points]();
    coord_values[1] = new double[num_points]();
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++) {
            coord_values[0][i * size + j] = 90.0  + 90.0 * j / size;
            coord_values[1][i * size + j] = -30.0 + 60.0 * i / size;
        } 
}
Grid_info_manager::~Grid_info_manager()
{
    delete coord_values[0];
    delete coord_values[1];
}
double** Grid_info_manager::get_grid_coord_values(int grid_id)
{
    return coord_values;
}
int Grid_info_manager::get_grid_num_points(int grid_id)
{
    return num_points;
}
void Grid_info_manager::get_grid_boundry(int grid_id, double* min_lat, double* max_lat, double* min_lon, double* max_lon)
{
    *min_lat = -30.0;
    *max_lat =  30.0;
    *min_lon =  90.0;
    *max_lon = 180.0;
}
int Grid_info_manager::get_polar_points(int grid_id, char polar)
{
    return 0;
}
bool Grid_info_manager::is_grid_cyclic(int grid_id)
{
    return false;
}
