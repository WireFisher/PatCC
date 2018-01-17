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
#include <tr1/unordered_map>
#include "merge_sort.h"
#include "ccpl_utils.h"

#define DEFAULT_EXPANGDING_RATIO (0.2)
#define PDLN_EXPECTED_EXPANDING_LOOP_TIMES (3)

#define PDLN_SPOLAR_MAX_LAT (-45.0)
#define PDLN_NPOLAR_MIN_LAT (45.0)

#define PDLN_MAX_ITER_COUNT (10)

#define PDLN_TOLERABLE_ERROR (0.0001)
#define PDLN_FLOAT_EQ_ERROR (1e-10)

#define PDLN_DECOMPOSE_COMMON_MODE (0)
#define PDLN_DECOMPOSE_SPOLAR_MODE (1)
#define PDLN_DECOMPOSE_NPOLAR_MODE (2)

#define PDLN_NODE_TYPE_COMMON PDLN_DECOMPOSE_COMMON_MODE
#define PDLN_NODE_TYPE_SPOLAR PDLN_DECOMPOSE_SPOLAR_MODE
#define PDLN_NODE_TYPE_NPOLAR PDLN_DECOMPOSE_NPOLAR_MODE

#define PDLN_DOUBLE_INVALID_VALUE ((double)0xDEADBEEFDEADBEEF)


//static inline bool is_cyclic(Boundry b)
//{
//    return std::abs(b.min_lon - 0.0) < PDLN_FLOAT_EQ_ERROR && std::abs(b.max_lon - 360.0) < PDLN_FLOAT_EQ_ERROR;
//}


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
    min_lat -= (max_lat - min_lat) * ratio * 0.5;
    max_lat += (max_lat - min_lat) * ratio * 0.5;
    min_lon -= (max_lon - min_lon) * ratio * 0.5;
    max_lon += (max_lon - min_lon) * ratio * 0.5;
    return *this;
}


void Boundry::legalize()
{
    min_lat = std::max(min_lat, -90.0);
    min_lon = std::max(min_lon, 0.0);
    max_lat = std::min(max_lat, 90.0);
    max_lon = std::min(max_lon, 360.0);
}


void Boundry::legalize(const Boundry *outer_boundry, bool is_cyclic)
{
    min_lat = std::max(min_lat, outer_boundry->min_lat);
    max_lat = std::min(max_lat, outer_boundry->max_lat);
    if(!is_cyclic) {
        min_lon = std::max(min_lon, outer_boundry->min_lon);
        max_lon = std::min(max_lon, outer_boundry->max_lon);
    }
}


Search_tree_node::Search_tree_node(Search_tree_node *parent, double *coord_value[2], int *global_index, int num_points, Boundry boundry) {
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
    this->local_cells_global_index = new int[num_points];
    memcpy(this->local_cells_coord[0], coord_value[0], num_points * sizeof(double));
    memcpy(this->local_cells_coord[1], coord_value[1], num_points * sizeof(double));
    memcpy(this->local_cells_global_index, global_index, num_points * sizeof(int));
    this->len_expanded_cells_coord_buf = 0;
    this->rotated_cells_coord[0] = this->rotated_cells_coord[1] = NULL;
    this->num_rotated_cells = 0;
    //this->expanded_cells_coord[0] = this->expanded_cells_coord[1] = NULL;

    this->num_local_kernel_cells = num_points;
    this->num_local_expanded_cells = 0;
    this->node_type = PDLN_NODE_TYPE_COMMON;
    this->triangulation = NULL;
}


Search_tree_node::~Search_tree_node()
{
    delete kernel_boundry;
    delete expanded_boundry;
    delete[] local_cells_coord[0];
    delete[] local_cells_coord[1];
    delete[] local_cells_global_index;
    if(rotated_kernel_boundry != NULL)
        delete rotated_kernel_boundry;
    if(rotated_expanded_boundry != NULL)
        delete rotated_expanded_boundry;
    //delete[] this->local_cells_global_index;
    if(first_child != NULL)
        delete first_child;
    if(second_child != NULL)
        delete second_child;
    if(third_child != NULL)
        delete third_child;
    if(triangulation != NULL)
        delete triangulation;
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


void Search_tree_node::generate_local_triangulation(bool is_cyclic)
{
    if(triangulation != NULL)
        delete triangulation;

    if(rotated_expanded_boundry != NULL) {
        triangulation = new Delaunay_Voronoi(num_local_kernel_cells + num_local_expanded_cells,
                                             rotated_cells_coord[PDLN_LON], rotated_cells_coord[PDLN_LAT], local_cells_global_index, false,
                                             rotated_expanded_boundry->min_lon, rotated_expanded_boundry->max_lon,
                                             rotated_expanded_boundry->min_lat, rotated_expanded_boundry->max_lat, NULL);

        double min_lon = 361.0, max_lon = -1.0;
        for (int i = 0; i < num_local_kernel_cells + num_local_expanded_cells; i++) {
            if(local_cells_coord[PDLN_LON][i] < min_lon) min_lon = local_cells_coord[PDLN_LON][i];
            if(local_cells_coord[PDLN_LON][i] > max_lon) max_lon = local_cells_coord[PDLN_LON][i];
        }
        double min_lat = 91.0, max_lat = -91.0;
        for (int i = 0; i < num_local_kernel_cells + num_local_expanded_cells; i++) {
            if(local_cells_coord[PDLN_LAT][i] < min_lat) min_lat = local_cells_coord[PDLN_LAT][i];
            if(local_cells_coord[PDLN_LAT][i] > max_lat) max_lat = local_cells_coord[PDLN_LAT][i];
        }

        double lon, head_lon, head_lat, tail_lon, tail_lat;
        lon = (max_lon + min_lon + 360.0) * 0.5;
        if (lon > 360.0) lon -= 360.0;
        if (std::abs(expanded_boundry->max_lat - 90.0) < PDLN_FLOAT_EQ_ERROR) {
            rotate_sphere_coordinate(lon, max_lat-0.1, head_lon, head_lat);
            rotate_sphere_coordinate(lon, min_lat,     tail_lon, tail_lat);
        }
        else {
            rotate_sphere_coordinate(lon, min_lat+0.1, head_lon, head_lat);
            rotate_sphere_coordinate(lon, max_lat,     tail_lon, tail_lat);
        }
        head_lon += 90.0;
        tail_lon += 90.0;
        if(head_lon > 360.0) head_lon -= 360.0;
        if(tail_lon > 360.0) tail_lon -= 360.0;

        //printf("lon: %lf max_lon: %lf, min_lon: %lf, max_lat: %lf, min_lat: %lf\n", lon, max_lon, min_lon, max_lat, min_lat);
        //printf("(%lf, %lf) -- (%lf, %lf)\n", head_lon, head_lat, tail_lon, tail_lat);
        std::vector<Triangle*> cyclic_triangles = triangulation->search_cyclic_triangles_for_rotated_grid(Point(head_lon, head_lat), Point(tail_lon, tail_lat));

        char filename[64];
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        snprintf(filename, 64, "log/cyclic_triangles_%d", rank);
        plot_triangles_info_file(filename, cyclic_triangles);

        triangulation->update_all_points_coord(local_cells_coord[PDLN_LON], local_cells_coord[PDLN_LAT], num_local_kernel_cells + num_local_expanded_cells);
        triangulation->correct_cyclic_triangles(cyclic_triangles, is_cyclic);
        triangulation->relegalize_all_triangles();
    }
    else {
        triangulation = new Delaunay_Voronoi(num_local_kernel_cells + num_local_expanded_cells,
                                             local_cells_coord[PDLN_LON], local_cells_coord[PDLN_LAT], local_cells_global_index, false,
                                             expanded_boundry->min_lon, expanded_boundry->max_lon,
                                             expanded_boundry->min_lat, expanded_boundry->max_lat, NULL);

    }
}

void Search_tree_node::decompose_with_certain_line(Midline midline, double *child_cells_coord[4], int *child_cells_idx[2], int child_num_cells[2])
{
    child_num_cells[0] = 0;
    child_num_cells[1] = 0;
    for(int i = 0; i < this->num_local_kernel_cells; i ++) {
        if(local_cells_coord[midline.type][i] < midline.value) {
            child_cells_coord[midline.type][child_num_cells[0]] = local_cells_coord[midline.type][i];
            child_cells_coord[(midline.type+1)%2][child_num_cells[0]] = local_cells_coord[(midline.type+1)%2][i];
            child_cells_idx[0][child_num_cells[0]++] = local_cells_global_index[i];
        }
        else {
            child_cells_coord[2+midline.type][child_num_cells[1]] = local_cells_coord[midline.type][i];
            child_cells_coord[2+(midline.type+1)%2][child_num_cells[1]] = local_cells_coord[(midline.type+1)%2][i];
            child_cells_idx[1][child_num_cells[1]++] = local_cells_global_index[i];
        }
    }
}


void Search_tree_node::decompose_iteratively(double *workloads, double *child_cells_coord[4], int *child_cells_idx[2], int child_num_cells[2],
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
        this->decompose_with_certain_line(midline, child_cells_coord, child_cells_idx, child_num_cells);
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
            this->decompose_with_certain_line(midline, child_cells_coord, child_cells_idx, child_num_cells);
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


void Search_tree_node::add_expanded_points(double *coord_value[2], int *global_idx, int num_points)
{
    double *tmp_coord_value[2];
    int *tmp_idx;
    if(num_local_expanded_cells + num_points > len_expanded_cells_coord_buf) {
        len_expanded_cells_coord_buf = num_local_expanded_cells + num_points * 4 * PDLN_EXPECTED_EXPANDING_LOOP_TIMES;
        tmp_coord_value[0] = new double[num_local_kernel_cells + len_expanded_cells_coord_buf];
        tmp_coord_value[1] = new double[num_local_kernel_cells + len_expanded_cells_coord_buf];
        tmp_idx = new int[num_local_kernel_cells + len_expanded_cells_coord_buf];
        memcpy(tmp_coord_value[0], local_cells_coord[0], sizeof(double) * (num_local_kernel_cells + num_local_expanded_cells));
        memcpy(tmp_coord_value[1], local_cells_coord[1], sizeof(double) * (num_local_kernel_cells + num_local_expanded_cells));
        memcpy(tmp_idx, local_cells_global_index, sizeof(int) * (num_local_kernel_cells + num_local_expanded_cells));
        delete[] local_cells_coord[0];
        delete[] local_cells_coord[1];
        delete[] local_cells_global_index;
        local_cells_coord[0] = tmp_coord_value[0];
        local_cells_coord[1] = tmp_coord_value[1];
        local_cells_global_index = tmp_idx;
        if(rotated_cells_coord[0] != NULL) {
            tmp_coord_value[0] = new double[num_local_kernel_cells + len_expanded_cells_coord_buf];
            tmp_coord_value[1] = new double[num_local_kernel_cells + len_expanded_cells_coord_buf];
            memcpy(tmp_coord_value[0], rotated_cells_coord[0], sizeof(double) * (num_local_kernel_cells + num_local_expanded_cells));
            memcpy(tmp_coord_value[1], rotated_cells_coord[1], sizeof(double) * (num_local_kernel_cells + num_local_expanded_cells));
            delete[] rotated_cells_coord[0];
            delete[] rotated_cells_coord[1];
            rotated_cells_coord[0] = tmp_coord_value[0];
            rotated_cells_coord[1] = tmp_coord_value[1];
        }
    }
    memcpy(local_cells_coord[0] + num_local_kernel_cells + num_local_expanded_cells, coord_value[0], sizeof(double) * num_points);
    memcpy(local_cells_coord[1] + num_local_kernel_cells + num_local_expanded_cells, coord_value[1], sizeof(double) * num_points);
    memcpy(local_cells_global_index + num_local_kernel_cells + num_local_expanded_cells, global_idx, sizeof(int) * num_points);
    num_local_expanded_cells += num_points;
    assert(num_local_expanded_cells >= num_points);
}


bool operator == (pair<Search_tree_node*, bool> p1, Search_tree_node* p2)
{
    return p1.first == p2;
}


void Search_tree_node::add_neighbors(vector<Search_tree_node*> neighbors)
{
    for(unsigned int i = 0; i < neighbors.size(); i++)
        if(find(this->neighbors.begin(), this->neighbors.end(), neighbors[i]) == this->neighbors.end())
            this->neighbors.push_back(pair<Search_tree_node*, bool>(neighbors[i], false));
}


/*
bool Search_tree_node::check_expanded_triangle_consistency()
{
    if(this->triangulation == NULL)
        return false;
}
*/


void Search_tree_node::generate_rotated_grid()
{
    if(rotated_cells_coord[0] == NULL) { /* first time to generate rotated grid */
        rotated_cells_coord[0] = new double[num_local_kernel_cells + len_expanded_cells_coord_buf];
        rotated_cells_coord[1] = new double[num_local_kernel_cells + len_expanded_cells_coord_buf];

        for(int i = 0; i < num_local_kernel_cells + num_local_expanded_cells; i++) {
            rotate_sphere_coordinate(local_cells_coord[PDLN_LON][i], local_cells_coord[PDLN_LAT][i], rotated_cells_coord[PDLN_LON][i], rotated_cells_coord[PDLN_LAT][i]);
            //int rank;
            //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            //if(rank == 0)
            //printf("coord: %lf, %lf\n", rotated_cells_coord[PDLN_LON][i], rotated_cells_coord[PDLN_LAT][i]);
            rotated_cells_coord[PDLN_LON][i] += 90; /* shift to avoid cyclic situation */
            if(rotated_cells_coord[PDLN_LON][i] >= 360.0) rotated_cells_coord[PDLN_LON][i] -= 360.0;
        }

        num_rotated_cells = num_local_kernel_cells + num_local_expanded_cells;

    }
    else {
        for(int i = num_rotated_cells; i < num_local_kernel_cells + num_local_expanded_cells; i++) {
            rotate_sphere_coordinate(local_cells_coord[PDLN_LON][i], local_cells_coord[PDLN_LAT][i], rotated_cells_coord[PDLN_LON][i], rotated_cells_coord[PDLN_LAT][i]);
            //int rank;
            //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            //if(rank == 0)
            //printf("coord: %lf, %lf\n", rotated_cells_coord[PDLN_LON][i], rotated_cells_coord[PDLN_LAT][i]);
            rotated_cells_coord[PDLN_LON][i] += 90;
            if(rotated_cells_coord[PDLN_LON][i] >= 360.0) rotated_cells_coord[PDLN_LON][i] -= 360.0;
        }
        num_rotated_cells = num_local_kernel_cells + num_local_expanded_cells;
    }

    /* recalculate expanded boundary */
    double top = -91.0, bot = 91.0, left = 540.0, right = -180.0;
    for(int i = 0; i < num_rotated_cells; i++) { //TODO: i can be started from non-zero
        if (rotated_cells_coord[PDLN_LON][i] < left)  left = rotated_cells_coord[PDLN_LON][i];
        if (rotated_cells_coord[PDLN_LON][i] > right) right = rotated_cells_coord[PDLN_LON][i];
        if (rotated_cells_coord[PDLN_LAT][i] < bot) bot = rotated_cells_coord[PDLN_LAT][i];
        if (rotated_cells_coord[PDLN_LAT][i] > top) top = rotated_cells_coord[PDLN_LAT][i];
    }
    if(rotated_expanded_boundry != NULL) {
        assert(rotated_expanded_boundry->min_lon >= left && rotated_expanded_boundry->max_lon <= right && rotated_expanded_boundry->min_lat >= bot && rotated_expanded_boundry->max_lat <= top);
        delete rotated_expanded_boundry;
    }
    rotated_expanded_boundry = new Boundry(left, right, bot, top);
}


Delaunay_grid_decomposition::Delaunay_grid_decomposition(int grid_id, Processing_resource *proc_info, int min_num_points_per_chunk)
{
    double **coord_values;
    int *global_index;
    Boundry boundry;
    int num_points;

    this->original_grid = grid_id;
    this->processing_info = proc_info;
    assert(this->processing_info != NULL);

    this->min_num_points_per_chunk = min_num_points_per_chunk;
    
    coord_values = grid_info_mgr->get_grid_coord_values(grid_id);
    num_points = grid_info_mgr->get_grid_num_points(grid_id);

    boundry.min_lat = 91.0;
    boundry.max_lat = -91.0;
    boundry.min_lon = 361.0;
    boundry.max_lon = -1.0;
    for(int i = 0; i < num_points; i++) {
        if(coord_values[PDLN_LON][i] < boundry.min_lon) boundry.min_lon = coord_values[PDLN_LON][i];
        if(coord_values[PDLN_LON][i] > boundry.max_lon) boundry.max_lon = coord_values[PDLN_LON][i];
    }
    for(int i = 0; i < num_points; i++) {
        if(coord_values[PDLN_LAT][i] < boundry.min_lat) boundry.min_lat = coord_values[PDLN_LAT][i];
        if(coord_values[PDLN_LAT][i] > boundry.max_lat) boundry.max_lat = coord_values[PDLN_LAT][i];
    }
    boundry.max_lon += 0.01;

    this->is_cyclic = grid_info_mgr->is_grid_cyclic(grid_id);
    //grid_info_mgr->get_grid_boundry(grid_id, &boundry.min_lat, &boundry.max_lat, &boundry.min_lon, &boundry.max_lon);
    this->processing_info->get_num_total_processing_units();

    global_index = new int[num_points];
    for(int i = 0; i < num_points; i++)
        global_index[i] = i;
    //for(int i = 0; i < num_points; i++)
    //    printf("(%lf, %lf)\n", coord_values[0][i], coord_values[1][i]);
    this->search_tree_root = new Search_tree_node(NULL, coord_values, global_index, num_points, boundry);
    this->active_processing_common_id = NULL;
    this->workloads = NULL;
    //this->initialze_workload();
    //
    delete global_index;
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

    if(ids.size() == 1) {
        this->workloads[ids[0]] = total_workload;
        return;
    }

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
    int *child_cells_index[2];
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
    child_cells_index[0] = new int[node->num_local_kernel_cells];
    child_cells_index[1] = new int[node->num_local_kernel_cells];
    
    node->decompose_iteratively(this->workloads, child_cells_coord, child_cells_index, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);

    node->first_child = new Search_tree_node(node, child_cells_coord,   child_cells_index[0], child_num_cells[0], child_boundry[0]);
    node->third_child = new Search_tree_node(node, child_cells_coord+2, child_cells_index[1], child_num_cells[1], child_boundry[1]);

    /* child_proc_id[0] can be modified by this->update_workloads */
    this->update_workloads(child_num_cells[0], child_proc_id[0]);
    this->update_workloads(child_num_cells[1], child_proc_id[1]);
    node->first_child->update_processing_units_id(child_proc_id[0]);
    node->third_child->update_processing_units_id(child_proc_id[1]);

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
    delete[] child_cells_index[0];
    delete[] child_cells_index[1];
    //TODO: optimize new delete
    
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("[Rank%d]x[ST-INFO-PRE] p: %p, first: %p, third: %p\n", rank, node, node->first_child, node->third_child);

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


void Delaunay_grid_decomposition::compute_common_boundry(Search_tree_node *tree_node1, Search_tree_node *tree_node2,
                                                         Point *boundry_head, Point *boundry_tail, Point *boundry2_head, Point *boundry2_tail)
{
    double coord_value[2][2];

    coord_value[0][PDLN_LAT] = coord_value[0][PDLN_LON] = coord_value[1][PDLN_LAT] = coord_value[1][PDLN_LON] = PDLN_DOUBLE_INVALID_VALUE;
    if(fabs(tree_node1->kernel_boundry->min_lat - tree_node2->kernel_boundry->max_lat) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon) < 
           std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon)) {
            coord_value[0][PDLN_LAT] = coord_value[1][PDLN_LAT] = tree_node2->kernel_boundry->max_lat;
            coord_value[0][PDLN_LON] = std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon);
            coord_value[1][PDLN_LON] = std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon);
        }
    }
    else if(fabs(tree_node1->kernel_boundry->max_lat - tree_node2->kernel_boundry->min_lat) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon) <
           std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon)) {
            coord_value[0][PDLN_LAT] = coord_value[1][PDLN_LAT] = tree_node1->kernel_boundry->max_lat;
            coord_value[0][PDLN_LON] = std::max(tree_node1->kernel_boundry->min_lon, tree_node2->kernel_boundry->min_lon);
            coord_value[1][PDLN_LON] = std::min(tree_node1->kernel_boundry->max_lon, tree_node2->kernel_boundry->max_lon);
        }
    }
    else if(fabs(tree_node1->kernel_boundry->min_lon - tree_node2->kernel_boundry->max_lon) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) <
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node2->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    else if(fabs(tree_node1->kernel_boundry->max_lon - tree_node2->kernel_boundry->min_lon) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) <
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node1->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    *boundry_head = Point(coord_value[0][PDLN_LON], coord_value[0][PDLN_LAT]);
    *boundry_tail = Point(coord_value[1][PDLN_LON], coord_value[1][PDLN_LAT]);

    coord_value[0][PDLN_LAT] = coord_value[0][PDLN_LON] = coord_value[1][PDLN_LAT] = coord_value[1][PDLN_LON] = PDLN_DOUBLE_INVALID_VALUE;
    if(fabs(fabs(tree_node1->kernel_boundry->min_lon - tree_node2->kernel_boundry->max_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) < 
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node2->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    else if(fabs(fabs(tree_node1->kernel_boundry->max_lon - tree_node2->kernel_boundry->min_lon) - 360.0) < PDLN_FLOAT_EQ_ERROR) {
        if(std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat) <
           std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat)) {
            coord_value[0][PDLN_LON] = coord_value[1][PDLN_LON] = tree_node1->kernel_boundry->max_lon;
            coord_value[0][PDLN_LAT] = std::max(tree_node1->kernel_boundry->min_lat, tree_node2->kernel_boundry->min_lat);
            coord_value[1][PDLN_LAT] = std::min(tree_node1->kernel_boundry->max_lat, tree_node2->kernel_boundry->max_lat);
        }
    }
    *boundry2_head = Point(coord_value[0][PDLN_LON], coord_value[0][PDLN_LAT]);
    *boundry2_tail = Point(coord_value[1][PDLN_LON], coord_value[1][PDLN_LAT]);
}


/* non-block */
void Delaunay_grid_decomposition::send_triangles_to_remote(int src_common_id, int dst_common_id, Triangle_Transport *triangles_buf, int num_triangles, int tag)
{
    MPI_Request request;
    if(processing_info->get_local_process_id() == processing_info->get_processing_unit(dst_common_id)->process_id)
        processing_info->send_to_local_thread(triangles_buf, num_triangles, sizeof(Triangle_Transport),
                                              processing_info->get_processing_unit(src_common_id)->thread_id,
                                              processing_info->get_processing_unit(dst_common_id)->thread_id, tag);
    else {
        MPI_Isend(triangles_buf, num_triangles*sizeof(Triangle_Transport), MPI_CHAR,
                  processing_info->get_processing_unit(dst_common_id)->process_id, 
                  tag, processing_info->get_mpi_comm(), &request);
        /*
        char filename[64];
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        snprintf(filename, 64, "log/sending_triangle_%dto%d", rank, processing_info->get_processing_unit(dst_common_id)->process_id);
        plot_triangles_info_file(filename, triangles_buf, num_triangles);
        */
    }
}


/* block */
int Delaunay_grid_decomposition::recv_triangles_from_remote(int src_common_id, int dst_common_id, Triangle_Transport *triangles_buf, int num_max_triangles, int tag)
{
    if(processing_info->get_local_process_id() == processing_info->get_processing_unit(src_common_id)->process_id)
        return processing_info->recv_from_local_thread(triangles_buf, num_max_triangles, sizeof(Triangle_Transport),
                                                       processing_info->get_processing_unit(src_common_id)->thread_id,
                                                       processing_info->get_processing_unit(dst_common_id)->thread_id, tag) / sizeof(Triangle_Transport);
    else {
        MPI_Status status;
        int count;
        MPI_Recv(triangles_buf, num_max_triangles*sizeof(Triangle_Transport), MPI_CHAR,
                 processing_info->get_processing_unit(src_common_id)->process_id, 
                 tag, processing_info->get_mpi_comm(), &status);
        MPI_Get_count(&status, MPI_CHAR, &count);
        assert(count%sizeof(Triangle_Transport) == 0);
        return count/sizeof(Triangle_Transport);
    }
}


namespace std   
{
    namespace tr1
    {
        template <>
        struct hash<Point>
        {
            std::size_t operator()(const Point &p) const
            {
                return (hash<double>()(p.x) ^ hash<double>()(p.y) << 1) >> 1;
            }
        };

        template <>  
        struct hash<Triangle_Transport>  
        {  
            std::size_t operator()(const Triangle_Transport &t) const  
            {  
                return hash<Point>()(t.v[0]) ^ hash<Point>()(t.v[1]) ^ hash<Point>()(t.v[2]);
            }
        }; 
    }
} 


bool Delaunay_grid_decomposition::check_triangles_consistency(Triangle_Transport *triangles1, Triangle_Transport *triangles2, int num_triangles)
{
    std::tr1::unordered_map<Triangle_Transport, bool> hash_table;

    if(num_triangles == 0)
        return true;

    for(int i = 0; i < num_triangles; i++) {
        if(hash_table.find(triangles1[i]) != hash_table.end())
            assert(false); //triangles1 has redundant triangles
        hash_table[triangles1[i]] = true;
    }

    for(int i = 0; i < num_triangles; i++) {
        if(hash_table.find(triangles2[i]) == hash_table.end())
            return false;
    }
    return true;
}


#define PDLN_COMM_TAG_MASK 0x0100
#define PDLN_SET_MASK(tag)       ((PDLN_COMM_TAG_MASK|tag)<<1|0)
#define PDLN_SET_MASK_EXTRA(tag) ((PDLN_COMM_TAG_MASK|tag)<<1|1)
bool Delaunay_grid_decomposition::check_leaf_node_triangulation_consistency(Search_tree_node *leaf_node, int iter)
{
    if(leaf_node->neighbors.size() == 0) {
        return false;
    }

    Triangle_Transport *local_triangle[32];
    Triangle_Transport *remote_triangle[32];
    Triangle_Transport *extra_local_triangle[32];
    Triangle_Transport *extra_remote_triangle[32];
    int triangle_buf_len;
    int num_local_triangle[32], num_remote_triangle[32];
    int num_extra_local_triangle[32], num_extra_remote_triangle[32];

    assert(leaf_node->processing_units_id.size() == 1);
    assert(leaf_node->neighbors.size() <= 32);
    assert(iter < PDLN_COMM_TAG_MASK);
    triangle_buf_len = search_tree_root->num_local_kernel_cells / processing_info->get_num_total_processing_units();
    
    for(unsigned int i = 0; i < leaf_node->neighbors.size(); i++) {
        local_triangle[i] = NULL;
        remote_triangle[i] = NULL;
        extra_local_triangle[i] = NULL;
        extra_remote_triangle[i] = NULL;
        num_local_triangle[i] = 0;
        num_remote_triangle[i] = 0;
        num_extra_local_triangle[i] = 0;
        num_extra_remote_triangle[i] = 0;
    }

    vector<MPI_Request*> waiting_list;

    /* send local triangles to neighbors */
    for(unsigned int i = 0; i < leaf_node->neighbors.size(); i++) {
        if(leaf_node->neighbors[i].second)
            continue;
        /* compute shared boundry of leaf_node and its neighbor */
        Point boundry_head, boundry_tail, boundry2_head, boundry2_tail;
        assert(leaf_node->neighbors[i].first->processing_units_id.size() == 1);
        compute_common_boundry(leaf_node, leaf_node->neighbors[i].first, &boundry_head, &boundry_tail, &boundry2_head, &boundry2_tail);
        /*
        if(leaf_node->processing_units_id[0] == 0)
        printf("[common boundary] %d -> %d: (%lf, %lf) -> (%lf, %lf)\n", leaf_node->processing_units_id[0],
                                                                         leaf_node->neighbors[i].first->processing_units_id[0],
                                                                         boundry_head.x, boundry_head.y,
                                                                         boundry_tail.x, boundry_tail.y);*/
        
        /* get triangles ,which the common boundry pass through, from leaf_node and its neighbor */
        if(boundry_head.x != PDLN_DOUBLE_INVALID_VALUE) { // Uncompleted determination
            local_triangle[i] = new Triangle_Transport[triangle_buf_len];
            remote_triangle[i] = new Triangle_Transport[triangle_buf_len];
            leaf_node->triangulation->get_triangles_intersecting_with_segment(boundry_head, boundry_tail, local_triangle[i], &num_local_triangle[i], triangle_buf_len);
            //printf("[%d]x[Send---] %d -> %d, num: %d, tag: %d\n", iter, leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0], num_local_triangle[i], PDLN_SET_MASK(iter));
            //send_triangles_to_remote(leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0],
            //                         local_triangle[i], num_local_triangle[i], PDLN_SET_MASK(iter));
            waiting_list.push_back(new MPI_Request);
            MPI_Isend(local_triangle[i], num_local_triangle[i]*sizeof(Triangle_Transport), MPI_CHAR,
                      processing_info->get_processing_unit(leaf_node->neighbors[i].first->processing_units_id[0])->process_id, 
                      PDLN_SET_MASK(iter), processing_info->get_mpi_comm(), waiting_list.back());
        }

        if(boundry2_head.x != PDLN_DOUBLE_INVALID_VALUE) { // Uncompleted determination
            printf("[common extra boundary] %d -> %d: (%lf, %lf) -> (%lf, %lf)\n", leaf_node->processing_units_id[0],
                                                                             leaf_node->neighbors[i].first->processing_units_id[0],
                                                                             boundry2_head.x, boundry2_head.y,
                                                                             boundry2_tail.x, boundry2_tail.y);
            extra_local_triangle[i] = new Triangle_Transport[triangle_buf_len];
            extra_remote_triangle[i] = new Triangle_Transport[triangle_buf_len];
            leaf_node->triangulation->get_triangles_intersecting_with_segment(boundry2_head, boundry2_tail, extra_local_triangle[i],
                                                                              &num_extra_local_triangle[i], triangle_buf_len);
            //printf("[%d]x[Send-EX] %d -> %d, num: %d, tag: %d\n", iter, leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0], num_local_triangle[i], PDLN_SET_MASK(iter));
            //send_triangles_to_remote(leaf_node->processing_units_id[0], leaf_node->neighbors[i].first->processing_units_id[0],
            //                          extra_local_triangle[i], num_extra_local_triangle[i], PDLN_SET_MASK_EXTRA(iter));
            waiting_list.push_back(new MPI_Request);
            MPI_Isend(extra_local_triangle[i], num_extra_local_triangle[i]*sizeof(Triangle_Transport), MPI_CHAR,
                      processing_info->get_processing_unit(leaf_node->neighbors[i].first->processing_units_id[0])->process_id, 
                      PDLN_SET_MASK_EXTRA(iter), processing_info->get_mpi_comm(), waiting_list.back());
        }
    }

    /* recv triangles from neighbors */
    for(unsigned int i = 0; i < leaf_node->neighbors.size(); i++) {
        if(leaf_node->neighbors[i].second)
            continue;
        Point boundry_head, boundry_tail, boundry2_head, boundry2_tail;
        compute_common_boundry(leaf_node, leaf_node->neighbors[i].first, &boundry_head, &boundry_tail, &boundry2_head, &boundry2_tail);

        if(boundry_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            assert(remote_triangle[i] != NULL);
            //printf("[%d]x[Recv---] %d -> %d, num: %d, tag: %d\n", iter, leaf_node->neighbors[i].first->processing_units_id[0], leaf_node->processing_units_id[0], triangle_buf_len, PDLN_SET_MASK(iter));
            num_remote_triangle[i] = recv_triangles_from_remote(leaf_node->neighbors[i].first->processing_units_id[0], leaf_node->processing_units_id[0],
                                                                remote_triangle[i], triangle_buf_len, PDLN_SET_MASK(iter));
            //printf("[%d]x[Recv+++] %d -> %d, num: %d, tag: %d\n", iter, leaf_node->neighbors[i].first->processing_units_id[0], leaf_node->processing_units_id[0], triangle_buf_len, PDLN_SET_MASK(iter));
        }

        if(boundry2_head.x != PDLN_DOUBLE_INVALID_VALUE) {
            assert(extra_remote_triangle[i] != NULL);
            //printf("[%d]x[Recv-EX] %d -> %d, num: %d, tag: %d\n", iter, leaf_node->neighbors[i].first->processing_units_id[0], leaf_node->processing_units_id[0], triangle_buf_len, PDLN_SET_MASK(iter));
            num_extra_remote_triangle[i] = recv_triangles_from_remote(leaf_node->neighbors[i].first->processing_units_id[0], leaf_node->processing_units_id[0],
                                                                      extra_remote_triangle[i], triangle_buf_len, PDLN_SET_MASK_EXTRA(iter));
        }
    }

    MPI_Status status;
    for(unsigned int i = 0; i < waiting_list.size(); i++)
        MPI_Wait(waiting_list[i], &status);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* do comparision */
    bool check_passed = true;
    for (unsigned int i = 0; i < leaf_node->neighbors.size(); i++) {
        if (leaf_node->neighbors[i].second)
            continue;
        if (local_triangle[i] != NULL) {
            if (num_local_triangle[i] != num_remote_triangle[i]) {
                char filename[64];
                snprintf(filename, 64, "log/boundary_triangle_local%d", rank);
                plot_triangles_info_file(filename, local_triangle[i], num_local_triangle[i]);
                snprintf(filename, 64, "log/boundary_triangle_remot%d", rank);
                plot_triangles_info_file(filename, remote_triangle[i], num_remote_triangle[i]);
                
                printf("[%d] checking consistency %d vs %d: number fault %d, %d\n", rank, leaf_node->processing_units_id[0],
                                                                                    leaf_node->neighbors[i].first->processing_units_id[0],
                                                                                    num_local_triangle[i], num_remote_triangle[i]);
                check_passed = false;
                continue;
            }
            if (!check_triangles_consistency(local_triangle[i], remote_triangle[i], num_local_triangle[i])) {
                char filename[64];
                snprintf(filename, 64, "log/boundary_triangle_local%d", rank);
                plot_triangles_info_file(filename, local_triangle[i], num_local_triangle[i]);
                snprintf(filename, 64, "log/boundary_triangle_remot%d", rank);
                plot_triangles_info_file(filename, remote_triangle[i], num_remote_triangle[i]);

                printf("[%d] checking consistency %d vs %d: triangle fault\n", rank, leaf_node->processing_units_id[0],
                                                                               leaf_node->neighbors[i].first->processing_units_id[0]);
                check_passed = false;
                continue;
            }
        }

        if (extra_local_triangle[i] != NULL) {
            if (num_extra_local_triangle[i] != num_extra_remote_triangle[i]) {
                printf("[%d] checking extra consistency %d vs %d: number fault %d, %d\n", rank, leaf_node->processing_units_id[0],
                                                                                          leaf_node->neighbors[i].first->processing_units_id[0],
                                                                                          num_local_triangle[i], num_remote_triangle[i]);
                check_passed = false;
                continue;
            }
            if (!check_triangles_consistency(extra_local_triangle[i], extra_remote_triangle[i], num_extra_local_triangle[i])) {
                printf("[%d] checking extra consistency %d vs %d: triangle fault\n", rank, leaf_node->processing_units_id[0], 
                                                                                     leaf_node->neighbors[i].first->processing_units_id[0]);
                check_passed = false;
                continue;
            }
        }
        leaf_node->neighbors[i].second = true;
    }

    if(check_passed) {
        printf("[%d] checking consistency %d: Pass\n", rank, leaf_node->processing_units_id[0]);
        return true;
    }
    else {
        //printf("[%d] checking consistency %d: Fail\n", rank, leaf_node->processing_units_id[0]);
        return false;
    }
}


void Delaunay_grid_decomposition::search_leaf_nodes_overlapping_with_region_recursively(Search_tree_node *node, Boundry region, vector<Search_tree_node*> &leaf_nodes_found)
{
    double *child_cells_coord[4];
    int *child_cells_index[2];
    int child_num_cells[2];
    Boundry child_boundry[2];
    vector<int> child_proc_id[2];

    if(node->processing_units_id.size() <= 0){
        assert(false);
        //print_whole_search_tree_info();
        //printf("node: %p\n", node);
    }
    assert(node->processing_units_id.size() > 0);
    if(node->processing_units_id.size() == 1) {
        if(two_regions_overlap(region, *node->kernel_boundry) ||
           two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->kernel_boundry) ||
           two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->kernel_boundry))
            leaf_nodes_found.push_back(node);
        return;
    }

    if(node->first_child == NULL && node->second_child == NULL && node->third_child == NULL) {
        for(int i = 0; i < 4; i++)
            child_cells_coord[i] = new double[node->num_local_kernel_cells];
        child_cells_index[0] = new int[node->num_local_kernel_cells];
        child_cells_index[1] = new int[node->num_local_kernel_cells];
        
        node->decompose_iteratively(this->workloads, child_cells_coord, child_cells_index, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_COMMON_MODE);
        assert(child_proc_id[0].size() > 0);

        node->first_child = new Search_tree_node(node, child_cells_coord,   child_cells_index[0], child_num_cells[0], child_boundry[0]);
        node->third_child = new Search_tree_node(node, child_cells_coord+2, child_cells_index[1], child_num_cells[1], child_boundry[1]);

        /* child_proc_id[0] can be modified by this->update_workloads */
        this->update_workloads(child_num_cells[0], child_proc_id[0]);
        this->update_workloads(child_num_cells[1], child_proc_id[1]);
        node->first_child->update_processing_units_id(child_proc_id[0]);
        assert(child_proc_id[1].size() > 0);
        node->third_child->update_processing_units_id(child_proc_id[1]);

        assert(node->first_child->processing_units_id.size() > 0);
        assert(node->third_child->processing_units_id.size() > 0);

        for(int i = 0; i < 4; i++)
            delete[] child_cells_coord[i];
        delete[] child_cells_index[0];
        delete[] child_cells_index[1];
        //TODO: optimize new delete
    }

    if(two_regions_overlap(region, *node->first_child->kernel_boundry) ||
       two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->first_child->kernel_boundry) ||
       two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->first_child->kernel_boundry))
        this->search_leaf_nodes_overlapping_with_region_recursively(node->first_child, region, leaf_nodes_found);
    if(node->second_child != NULL)
        if(two_regions_overlap(region, *node->second_child->kernel_boundry) ||
           two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->second_child->kernel_boundry) ||
           two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->second_child->kernel_boundry))
            this->search_leaf_nodes_overlapping_with_region_recursively(node->second_child, region, leaf_nodes_found);
    if(two_regions_overlap(region, *node->third_child->kernel_boundry) ||
       two_regions_overlap(Boundry(region.min_lon + 360.0, region.max_lon + 360.0, region.min_lat, region.max_lat), *node->third_child->kernel_boundry) ||
       two_regions_overlap(Boundry(region.min_lon - 360.0, region.max_lon - 360.0, region.min_lat, region.max_lat), *node->third_child->kernel_boundry))
        this->search_leaf_nodes_overlapping_with_region_recursively(node->third_child, region, leaf_nodes_found);
}


int Delaunay_grid_decomposition::assign_polars(bool assign_south_polar, bool assign_north_polar)
{
    double *child_cells_coord[4];
    int *child_cells_index[2];
    Boundry child_boundry[2];
    int child_num_cells[2];
    vector<int> child_proc_id[2];
    Midline midline;
    
    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[this->search_tree_root->num_local_kernel_cells];
    child_cells_index[0] = new int[this->search_tree_root->num_local_kernel_cells];
    child_cells_index[1] = new int[this->search_tree_root->num_local_kernel_cells];

    if(assign_south_polar) {
        this->search_tree_root->decompose_iteratively(this->workloads, child_cells_coord, child_cells_index, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_SPOLAR_MODE);
        if(child_boundry[0].max_lat > PDLN_SPOLAR_MAX_LAT || this->search_tree_root->processing_units_id.size() == 1) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_SPOLAR_MAX_LAT;
            this->search_tree_root->decompose_with_certain_line(midline, child_cells_coord, child_cells_index, child_num_cells);

            if(child_num_cells[0] < this->min_num_points_per_chunk) {
                for(int i = 0; i < 4; i++)
                    delete[] child_cells_coord[i];
                delete[] child_cells_index[0];
                delete[] child_cells_index[1];
                return 1;
            }
            child_boundry[0].max_lat = PDLN_SPOLAR_MAX_LAT;
            child_boundry[1].min_lat = PDLN_SPOLAR_MAX_LAT;
            this->workloads[this->search_tree_root->processing_units_id[0]] -= child_num_cells[0];
            child_proc_id[1].insert(child_proc_id[1].begin(), this->search_tree_root->processing_units_id[0]);
        }
        this->search_tree_root->first_child  = new Search_tree_node(this->search_tree_root, child_cells_coord,   child_cells_index[0], child_num_cells[0], child_boundry[0]);
        this->search_tree_root->second_child = new Search_tree_node(this->search_tree_root, child_cells_coord+2, child_cells_index[1], child_num_cells[1], child_boundry[1]);

        this->search_tree_root->first_child->node_type = PDLN_NODE_TYPE_SPOLAR;

        //this->update_workloads(child_num_cells[0], child_proc_id[0]);
        this->update_workloads(child_num_cells[1], child_proc_id[1]);
        this->search_tree_root->first_child->update_processing_units_id(child_proc_id[0]);
        this->search_tree_root->second_child->update_processing_units_id(child_proc_id[1]);

        this->current_tree_node = this->search_tree_root->second_child;
        
        if(this->have_local_processing_units_id(this->search_tree_root->first_child->processing_units_id))
            this->local_leaf_nodes.push_back(this->search_tree_root->first_child);
    }
    
    if(assign_north_polar) {
        this->current_tree_node->decompose_iteratively(this->workloads, child_cells_coord, child_cells_index, child_num_cells, child_boundry, child_proc_id, PDLN_DECOMPOSE_NPOLAR_MODE);
        if(child_boundry[1].min_lat < PDLN_NPOLAR_MIN_LAT || this->current_tree_node->processing_units_id.size() == 1) {
            midline.type = PDLN_LAT;
            midline.value = PDLN_NPOLAR_MIN_LAT;
            this->current_tree_node->decompose_with_certain_line(midline, child_cells_coord, child_cells_index, child_num_cells);

            if(child_num_cells[1] < this->min_num_points_per_chunk) {
                for(int i = 0; i < 4; i++)
                    delete[] child_cells_coord[i];
                delete[] child_cells_index[0];
                delete[] child_cells_index[1];
                return 1;
            }
            child_boundry[0].max_lat = PDLN_NPOLAR_MIN_LAT;
            child_boundry[1].min_lat = PDLN_NPOLAR_MIN_LAT;
            this->workloads[this->search_tree_root->processing_units_id.back()] -= child_num_cells[1];
            child_proc_id[0].push_back(this->search_tree_root->processing_units_id.back());
        }
        if(this->search_tree_root->second_child != NULL)
            delete this->search_tree_root->second_child;
        this->search_tree_root->second_child = new Search_tree_node(this->search_tree_root, child_cells_coord,   child_cells_index[0], child_num_cells[0], child_boundry[0]);
        this->search_tree_root->third_child  = new Search_tree_node(this->search_tree_root, child_cells_coord+2, child_cells_index[1], child_num_cells[1], child_boundry[1]);

        this->search_tree_root->third_child->node_type = PDLN_NODE_TYPE_NPOLAR;

        this->update_workloads(child_num_cells[0], child_proc_id[0]);
        //this->update_workloads(child_num_cells[1], child_proc_id[1]);
        this->search_tree_root->second_child->update_processing_units_id(child_proc_id[0]);
        this->search_tree_root->third_child->update_processing_units_id(child_proc_id[1]);

        this->current_tree_node = this->search_tree_root->second_child;

        if(this->have_local_processing_units_id(this->search_tree_root->third_child->processing_units_id))
            this->local_leaf_nodes.push_back(this->search_tree_root->third_child);
    }

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
    delete[] child_cells_index[0];
    delete[] child_cells_index[1];
    return 0;
}


void Delaunay_grid_decomposition::assign_cyclic_grid_for_single_unit()
{
    double *child_cells_coord[4];
    int *child_cells_index[2];
    Boundry child_boundry[2];
    int child_num_cells[2];
    vector<int> child_proc_id[2];
    Midline midline;

    for(int i = 0; i < 4; i++)
        child_cells_coord[i] = new double[this->current_tree_node->num_local_kernel_cells];
    child_cells_index[0] = new int[this->search_tree_root->num_local_kernel_cells];
    child_cells_index[1] = new int[this->search_tree_root->num_local_kernel_cells];

    assert(this->current_tree_node->processing_units_id.size() == 1);
    midline.type = PDLN_LON;
    midline.value = 180.0;
    child_boundry[0] = child_boundry[1] = *this->current_tree_node->kernel_boundry;
    child_boundry[0].max_lon = child_boundry[1].min_lon = 180.0;
    child_proc_id[0] = child_proc_id[1] = this->current_tree_node->processing_units_id;
    this->current_tree_node->decompose_with_certain_line(midline, child_cells_coord, child_cells_index, child_num_cells);
    this->workloads[this->current_tree_node->processing_units_id[0]] -= child_num_cells[0] + child_num_cells[1];

    this->current_tree_node->first_child = new Search_tree_node(this->current_tree_node, child_cells_coord,   child_cells_index[0], child_num_cells[0], child_boundry[0]);
    this->current_tree_node->third_child = new Search_tree_node(this->current_tree_node, child_cells_coord+2, child_cells_index[1], child_num_cells[1], child_boundry[1]);
    this->current_tree_node->first_child->update_processing_units_id(child_proc_id[0]);
    this->current_tree_node->third_child->update_processing_units_id(child_proc_id[1]);

    assert(this->have_local_processing_units_id(child_proc_id[0]));
    this->local_leaf_nodes.push_back(this->current_tree_node->first_child);
    this->local_leaf_nodes.push_back(this->current_tree_node->third_child);

    for(int i = 0; i < 4; i++)
        delete[] child_cells_coord[i];
    delete[] child_cells_index[0];
    delete[] child_cells_index[1];
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
    if(this->assign_polars(std::abs(search_tree_root->kernel_boundry->min_lat - -90.0) < PDLN_FLOAT_EQ_ERROR && num_south_polar < 2,
                           std::abs(search_tree_root->kernel_boundry->max_lat -  90.0) < PDLN_FLOAT_EQ_ERROR && num_north_polar < 2))
        return 1;

    if(this->current_tree_node->processing_units_id.size() == 1 && is_cyclic) {
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


vector<Search_tree_node*> Delaunay_grid_decomposition::search_points_in_region(Boundry region, double *coord_values[2], int *global_idx, int *num_points_found)
{
    vector<Search_tree_node*> leaf_nodes_found;
    
    *num_points_found = 0;
    if(fabs(region.min_lat - region.max_lat) < PDLN_FLOAT_EQ_ERROR || fabs(region.min_lon - region.max_lon) < PDLN_FLOAT_EQ_ERROR)
        return leaf_nodes_found;

    search_leaf_nodes_overlapping_with_region_recursively(this->search_tree_root, region, leaf_nodes_found);

    double left_min_lon = region.min_lon - 360.0;
    double left_max_lon = region.max_lon - 360.0;
    double right_min_lon = region.min_lon + 360.0;
    double right_max_lon = region.max_lon + 360.0;

    //TODO: optimize if out of loop
    if((region.max_lon <= 360.0 && region.max_lon >= 0.0) || (region.min_lon <= 360.0 && region.min_lon >= 0.0))
        for(unsigned i = 0; i < leaf_nodes_found.size(); i++)
            for(int j = 0; j < leaf_nodes_found[i]->num_local_kernel_cells; j++)
                if (leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] < region.max_lon && leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] >= region.min_lon &&
                    leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] < region.max_lat && leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] >= region.min_lat) {
                    coord_values[PDLN_LON][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j];
                    coord_values[PDLN_LAT][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j];
                    global_idx[(*num_points_found)++] = leaf_nodes_found[i]->local_cells_global_index[j];
                }

    if((left_max_lon <= 360.0 && left_max_lon >= 0.0) || (left_min_lon <= 360.0 && left_min_lon >= 0.0))
        for(unsigned i = 0; i < leaf_nodes_found.size(); i++)
            for(int j = 0; j < leaf_nodes_found[i]->num_local_kernel_cells; j++)
                if (leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] < left_max_lon   && leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] >= left_min_lon &&
                    leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] < region.max_lat && leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] >= region.min_lat) {
                    coord_values[PDLN_LON][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] + 360.0;
                    coord_values[PDLN_LAT][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j];
                    global_idx[(*num_points_found)++] = leaf_nodes_found[i]->local_cells_global_index[j];
                }

    if((right_max_lon <= 360.0 && right_max_lon >= 0.0) || (right_min_lon <= 360.0 && right_min_lon >= 0.0))
        for(unsigned i = 0; i < leaf_nodes_found.size(); i++)
            for(int j = 0; j < leaf_nodes_found[i]->num_local_kernel_cells; j++)
                if (leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] < right_max_lon  && leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] >= right_min_lon &&
                    leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] < region.max_lat && leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j] >= region.min_lat) {
                    coord_values[PDLN_LON][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LON][j] - 360.0;
                    coord_values[PDLN_LAT][*num_points_found] = leaf_nodes_found[i]->local_cells_coord[PDLN_LAT][j];
                    global_idx[(*num_points_found)++] = leaf_nodes_found[i]->local_cells_global_index[j];
                }

    return leaf_nodes_found;
}


int Delaunay_grid_decomposition::expand_tree_node_boundry(Search_tree_node* tree_node, double expanding_ratio)
{
    Boundry old_boundry = *tree_node->expanded_boundry;
    Boundry sub_rectangles[4];
    double *expanded_cells_coord[2];
    int *expanded_index;
    int num_points_found;

    if(tree_node->node_type == PDLN_NODE_TYPE_COMMON)
        *tree_node->expanded_boundry = *tree_node->expanded_boundry * expanding_ratio;
    else if(tree_node->node_type == PDLN_NODE_TYPE_SPOLAR)
        tree_node->expanded_boundry->max_lat += (tree_node->expanded_boundry->max_lat - tree_node->expanded_boundry->min_lat) * expanding_ratio * 2;
    else if(tree_node->node_type == PDLN_NODE_TYPE_NPOLAR)
        tree_node->expanded_boundry->min_lat -= (tree_node->expanded_boundry->max_lat - tree_node->expanded_boundry->min_lat) * expanding_ratio * 2;

    tree_node->expanded_boundry->legalize(search_tree_root->kernel_boundry, is_cyclic);
    //if(tree_node->processing_units_id[0] == 0)
    //    printf("old: %lf, expanded: %lf\n", old_boundry.max_lat, tree_node->expanded_boundry->max_lat);

    /*
    if(tree_node->node_type == PDLN_NODE_TYPE_COMMON && old_boundry.max_lat >= 90 && tree_node->expanded_boundry->max_lat >= 90)
        return 1;
    if(tree_node->node_type == PDLN_NODE_TYPE_COMMON && old_boundry.min_lat <= -90 && tree_node->expanded_boundry->min_lat <= -90)
        return 1;
    */

    transform_into_rectangle(old_boundry, *tree_node->expanded_boundry, sub_rectangles);

    expanded_cells_coord[0] = new double[search_tree_root->num_local_kernel_cells]; //FIXME: buf too large
    expanded_cells_coord[1] = new double[search_tree_root->num_local_kernel_cells];
    expanded_index = new int[search_tree_root->num_local_kernel_cells];
    for(int i = 0; i < 4; i++){
        num_points_found = 0;
        vector<Search_tree_node*> leaf_nodes_found;
        leaf_nodes_found = search_points_in_region(sub_rectangles[i], expanded_cells_coord, expanded_index, &num_points_found);
        tree_node->add_expanded_points(expanded_cells_coord, expanded_index, num_points_found);
        tree_node->add_neighbors(leaf_nodes_found);
    }

    delete[] expanded_cells_coord[0];
    delete[] expanded_cells_coord[1];
    delete[] expanded_index;
    return 0;
}


bool Search_tree_node::check_if_all_outer_edge_out_of_kernel_boundry(Boundry *origin_grid_boundry, bool is_cyclic)
{
    double left, right, top, bot;
    if(triangulation == NULL)
        return false;

    top = fabs(kernel_boundry->max_lat - origin_grid_boundry->max_lat) < PDLN_FLOAT_EQ_ERROR ?
          (kernel_boundry->max_lat+kernel_boundry->min_lat)/2.0 : kernel_boundry->max_lat; //Make the value small enough to pass the check. This is not strict. May need FIXME.
    bot = fabs(kernel_boundry->min_lat - origin_grid_boundry->min_lat) < PDLN_FLOAT_EQ_ERROR ?
          (kernel_boundry->max_lat+kernel_boundry->min_lat)/2.0 : kernel_boundry->min_lat;
    if(is_cyclic) {
        left = kernel_boundry->max_lon;
        right = kernel_boundry->min_lon;
    }
    else {
        left = fabs(kernel_boundry->max_lon - origin_grid_boundry->max_lon) < PDLN_FLOAT_EQ_ERROR ?
               (kernel_boundry->max_lon+kernel_boundry->min_lon)/2.0 : kernel_boundry->max_lon;
        right = fabs(kernel_boundry->min_lon - origin_grid_boundry->min_lon) < PDLN_FLOAT_EQ_ERROR ?
                (kernel_boundry->max_lon+kernel_boundry->min_lon)/2.0 : kernel_boundry->min_lon;
    } 
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("[%d] checking boundary: %lf, %lf, %lf, %lf\n", rank, top, bot, left, right);

    //return triangulation->check_if_all_outer_edge_out_of_region(left, right, bot, top);
    if(triangulation->check_if_all_outer_edge_out_of_region(left, right, bot, top))
        return true;
    else {
        printf("checkfailed for %d\n", processing_units_id[0]);
        return false;
    }
}


/* 
 * Return Values: 0 - Success
 *                1 - Failed in expanding
 *                2 - Fail, polar decomp's expanded_boundry exceeded threshold
 */
int Delaunay_grid_decomposition::generate_trianglulation_for_local_decomp()
{
    //TODO: openmp parallel
    bool expand_fail = false;

    if(local_leaf_nodes.size() > 1)
        printf("local_leaf_nodes.size(): %lu\n", local_leaf_nodes.size());
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        int ret;
        double expanding_ratio = DEFAULT_EXPANGDING_RATIO;
        int iter = 0;
        while(!check_leaf_node_triangulation_consistency(local_leaf_nodes[i], iter) ||
              !local_leaf_nodes[i]->check_if_all_outer_edge_out_of_kernel_boundry(search_tree_root->kernel_boundry, is_cyclic)) {
            /* Actually, we don't have confidence in making the loop iterate only once by tuning. TODO: try to delete "don't" */
            ret = expand_tree_node_boundry(local_leaf_nodes[i], expanding_ratio);
            /*
            printf("[%d]ID:%d expanded boundary (%lf, %lf, %lf, %lf)\n", iter, local_leaf_nodes[0]->processing_units_id[0], local_leaf_nodes[0]->expanded_boundry->min_lon,
                                                                       local_leaf_nodes[0]->expanded_boundry->max_lon, local_leaf_nodes[0]->expanded_boundry->min_lat,
                                                                       local_leaf_nodes[0]->expanded_boundry->max_lat);
            
            printf("[%d]x[INFO] ID: %d, (%lf, %lf, %lf, %lf), Neighbors: [", iter, local_leaf_nodes[0]->processing_units_id[0], local_leaf_nodes[0]->kernel_boundry->min_lon,
                                                                       local_leaf_nodes[0]->kernel_boundry->max_lon, local_leaf_nodes[0]->kernel_boundry->min_lat,
                                                                       local_leaf_nodes[0]->kernel_boundry->max_lat);
            for(unsigned int j = 0; j < local_leaf_nodes[0]->neighbors.size(); j++)
                printf("%d=%p, ", local_leaf_nodes[0]->neighbors[j].first->processing_units_id[0], local_leaf_nodes[0]->neighbors[j].first);
            printf("]\n");
            
            printf("===========%d===========\n", iter);
            //print_whole_search_tree_info();*/

            int rets[512];
            assert(processing_info->get_num_total_processes() <= 512);
            memset(rets, 0, 512*sizeof(int));
            MPI_Allgather(&ret, 1, MPI_INT, rets, 1, MPI_INT, processing_info->get_mpi_comm());
            for (unsigned j = 0; j < 512; j++)
                if (rets[j] != 0) {
                    expand_fail = true;
                    break;
                }

            if(expand_fail)
                break;

            if (local_leaf_nodes[i]->node_type != PDLN_NODE_TYPE_COMMON) {
                bool spolar_success = true, npolar_success = true;

                if (local_leaf_nodes[i]->node_type == PDLN_NODE_TYPE_SPOLAR && local_leaf_nodes[i]->expanded_boundry->max_lat > 0.0)
                    spolar_success = false;
                if (local_leaf_nodes[i]->node_type == PDLN_NODE_TYPE_NPOLAR && local_leaf_nodes[i]->expanded_boundry->min_lat < 0.0)
                    npolar_success = false;

                if (spolar_success && npolar_success) {
                    /* if MPI_Ibcast is supported
                    MPI_Request req[2];
                    MPI_Status status;
                    MPI_Ibcast(&spolar_success, 1*sizeof(bool), MPI_CHAR,
                               processing_info->get_processing_unit(search_tree_root->first_child->processing_units_id[0])->process_id, 
                               processing_info->get_mpi_comm(), req);
                    MPI_Ibcast(&npolar_success, 1*sizeof(bool), MPI_CHAR,
                               processing_info->get_processing_unit(search_tree_root->third_child->processing_units_id[0])->process_id, 
                               processing_info->get_mpi_comm(), req+1);

                    local_leaf_nodes[i]->generate_rotated_grid();

                    MPI_Wait(req, &status);
                    MPI_Wait(req+1, &status);
                    if (!spolar_success || !npolar_success)
                        return 2;
                    */
                    MPI_Bcast(&spolar_success, 1*sizeof(bool), MPI_CHAR,
                               processing_info->get_processing_unit(search_tree_root->first_child->processing_units_id[0])->process_id, 
                               processing_info->get_mpi_comm());
                    MPI_Bcast(&npolar_success, 1*sizeof(bool), MPI_CHAR,
                               processing_info->get_processing_unit(search_tree_root->third_child->processing_units_id[0])->process_id, 
                               processing_info->get_mpi_comm());
                    if (!spolar_success || !npolar_success)
                        return 2;

                    local_leaf_nodes[i]->generate_rotated_grid();
                }
                else {
                    MPI_Bcast(&spolar_success, 1*sizeof(bool), MPI_CHAR,
                              processing_info->get_processing_unit(search_tree_root->first_child->processing_units_id[0])->process_id, 
                              processing_info->get_mpi_comm());
                    MPI_Bcast(&npolar_success, 1*sizeof(bool), MPI_CHAR,
                              processing_info->get_processing_unit(search_tree_root->third_child->processing_units_id[0])->process_id, 
                              processing_info->get_mpi_comm());
                    return 2;
                }
            }
            else {
                bool spolar_success = true, npolar_success = true;
                MPI_Bcast(&spolar_success, 1*sizeof(bool), MPI_CHAR,
                          processing_info->get_processing_unit(search_tree_root->first_child->processing_units_id[0])->process_id, 
                          processing_info->get_mpi_comm());
                MPI_Bcast(&npolar_success, 1*sizeof(bool), MPI_CHAR,
                          processing_info->get_processing_unit(search_tree_root->third_child->processing_units_id[0])->process_id, 
                          processing_info->get_mpi_comm());
                if(!spolar_success || !npolar_success)
                    return 2;
            }

            local_leaf_nodes[i]->generate_local_triangulation(is_cyclic);
            expanding_ratio = DEFAULT_EXPANGDING_RATIO + 0.1;
            iter++;
            //if(iter == 2) {
                //plot_local_triangles("log/chunk");
                //MPI_Barrier(MPI_COMM_WORLD);
                //exit(1);
                //break;
            //}
        }
    }
    if(expand_fail)
        return 1;
    return 0;
}

int Delaunay_grid_decomposition::generate_trianglulation_for_whole_grid()
{
    assert(false);
    return 0;
}


void Delaunay_grid_decomposition::plot_local_triangles(const char *perfix)
{
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        char filename[64];
        snprintf(filename, 64, "%s_%d.png", perfix, local_leaf_nodes[i]->processing_units_id[0]);
        local_leaf_nodes[i]->triangulation->plot_into_file(filename, local_leaf_nodes[i]->kernel_boundry->min_lon,
                                                                     local_leaf_nodes[i]->kernel_boundry->max_lon,
                                                                     local_leaf_nodes[i]->kernel_boundry->min_lat,
                                                                     local_leaf_nodes[i]->kernel_boundry->max_lat);
    }
}


void Delaunay_grid_decomposition::print_tree_node_info_recursively(Search_tree_node *node)
{
    if(node->processing_units_id.size() == 1){
        printf("[ID%d]x[ST-Info] LEAF %p\n", local_leaf_nodes[0]->processing_units_id[0], node);
        return;
    }

    printf("[ID%d]x[ST-Info] %p: %p, %p, %p\n", local_leaf_nodes[0]->processing_units_id[0], node, node->first_child, node->second_child, node->third_child);
    if(node->first_child)
        print_tree_node_info_recursively(node->first_child);
    if(node->second_child)
        print_tree_node_info_recursively(node->second_child);
    if(node->third_child)
        print_tree_node_info_recursively(node->third_child);
}


void Delaunay_grid_decomposition::print_whole_search_tree_info()
{
    printf("[ID%d]x[ST-Info] ROOT %p\n", local_leaf_nodes[0]->processing_units_id[0], search_tree_root);
    print_tree_node_info_recursively(search_tree_root);
}



static inline void swap(Point *p1, Point *p2)
{
    Point tmp = *p1;
    *p1 = *p2;
    *p2 = tmp;
}

int compare_v2(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[2].id < t2.v[2].id) return -1;
    if(t1.v[2].id > t2.v[2].id) return  1;
    return 0;
}

int compare_v1(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[1].id < t2.v[1].id) return -1;
    if(t1.v[1].id > t2.v[1].id) return  1;
    return 0;
}

int compare_v0(const void* a, const void* b)
{
    Triangle_Transport t1 = *(const Triangle_Transport*)a;
    Triangle_Transport t2 = *(const Triangle_Transport*)b;

    if(t1.v[0].id < t2.v[0].id) return -1;
    if(t1.v[0].id > t2.v[0].id) return  1;
    return 0;
}

int compare_int(const void* a, const void* b)
{
    int t1 = *(const int*)a;
    int t2 = *(const int*)b;

    if(t1 < t2) return -1;
    if(t1 > t2) return  1;
    return 0;
}

static void radix_sort(Triangle_Transport *triangles, int num_triangles)
{
    assert(sizeof(Triangle_Transport) > sizeof(void *)/2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_v2);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_v1);
    merge_sort(triangles, num_triangles, sizeof(Triangle_Transport), compare_v0);
}

void Delaunay_grid_decomposition::save_ordered_triangles_into_file(Triangle_Transport *triangles, int num_triangles)
{
    int i, j;
    for(i = 0; i < num_triangles; i++) {
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
        if(triangles[i].v[1].id > triangles[i].v[2].id) swap(&triangles[i].v[1], &triangles[i].v[2]);
        if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
    }

    radix_sort(triangles, num_triangles);

    for(i = 0, j = 1; j < num_triangles; j++) {
        if(triangles[i].v[0].id == triangles[j].v[0].id &&
           triangles[i].v[1].id == triangles[j].v[1].id &&
           triangles[i].v[2].id == triangles[j].v[2].id)
            continue;
        else
            triangles[++i] = triangles[j];
    }
    int num_different_triangles = i + 1;
    
    FILE *fp = fopen("log/global_triangles", "w");
    for(i = 0; i < num_different_triangles; i++)
        fprintf(fp, "%d, %d, %d\n", triangles[i].v[0].id, triangles[i].v[1].id, triangles[i].v[2].id);
    fclose(fp);

    //char filename[64];
    //snprintf(filename, 64, "log/image_global_triangles%d");
    //plot_triangles_info_file("log/image_global_triangles", triangles, num_different_triangles);

}


#define PDLN_MERGE_TAG_MASK 0x0200
void Delaunay_grid_decomposition::merge_all_triangles()
{
    int local_buf_len = 0;
    int num_triangles, num_local_triangles;
    Triangle_Transport *local_triangles;

    /* Let n be the number of points, if there are b vertices on the convex hull,
     * then any triangulation of the points has at most 2n − 2 − b triangles,
     * plus one exterior face */
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++)
        local_buf_len += local_leaf_nodes[i]->num_local_kernel_cells * 3; 

    local_triangles = new Triangle_Transport[local_buf_len];
    num_local_triangles = 0;
    for(unsigned int i = 0; i < local_leaf_nodes.size(); i++) {
        local_leaf_nodes[i]->triangulation->get_triangles_in_region(local_leaf_nodes[i]->kernel_boundry->min_lon, local_leaf_nodes[i]->kernel_boundry->max_lon,
                                                                    local_leaf_nodes[i]->kernel_boundry->min_lat, local_leaf_nodes[i]->kernel_boundry->max_lat,
                                                                    local_triangles + num_local_triangles, &num_triangles, local_buf_len - num_local_triangles); // This function call is un-class-elegant
        num_local_triangles += num_triangles;
    }


    MPI_Barrier(processing_info->get_mpi_comm());

    if(processing_info->get_local_process_id() == 0) {
        int *num_remote_triangles = new int[processing_info->get_num_total_processes()];
        int remote_buf_len = 0;
        Triangle_Transport *remote_triangles;
        MPI_Status status;

        for(int i = 1; i < processing_info->get_num_total_processes(); i++)
            MPI_Recv(&num_remote_triangles[i], 1, MPI_INT, i, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm(), &status);
        for(int i = 1; i < processing_info->get_num_total_processes(); i++) {
            assert(num_remote_triangles[i] > min_num_points_per_chunk/2);
            remote_buf_len += num_remote_triangles[i];
        }
        remote_triangles = new Triangle_Transport[remote_buf_len + num_local_triangles];

        int count = 0;
        for(int i = 1; i < processing_info->get_num_total_processes(); i++) {
            MPI_Recv(remote_triangles + count, num_remote_triangles[i] * sizeof(Triangle_Transport), MPI_CHAR, i, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm(), &status);
            int tmp_count;
            MPI_Get_count(&status, MPI_CHAR, &tmp_count);
            
            //char filename[64];
            //snprintf(filename, 64, "log/process_local_triangles%d", i);
            //plot_triangles_info_file(filename, remote_triangles+count, tmp_count/sizeof(Triangle_Transport));

            //for(int j = 0; j < tmp_count/sizeof(Triangle_Transport); j++)
            //    printf("%d, %d, %d\n", (remote_triangles + count)[j].v[0].id, (remote_triangles + count)[j].v[1].id, (remote_triangles + count)[j].v[2].id);
            //printf("==============\n");
#ifdef DEBUG
            assert(tmp_count % sizeof(Triangle_Transport) == 0);
#endif
            count += tmp_count / sizeof(Triangle_Transport);
        }
        assert(count == remote_buf_len);
        memcpy(remote_triangles + remote_buf_len, local_triangles, num_local_triangles * sizeof(Triangle_Transport));
        save_ordered_triangles_into_file(remote_triangles, remote_buf_len + num_local_triangles);
    }
    else {
        MPI_Send(&num_local_triangles, 1, MPI_INT, 0, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm());
        MPI_Send(local_triangles, num_local_triangles * sizeof(Triangle_Transport), MPI_CHAR, 0, PDLN_MERGE_TAG_MASK, processing_info->get_mpi_comm());
    }
}


static double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


Grid_info_manager::Grid_info_manager()
{
    int size = 300;

    num_points = size * size;
    coord_values[0] = new double[num_points]();
    coord_values[1] = new double[num_points]();
    srand(0);
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++) {
            coord_values[0][i * size + j] =  0.0  + 360.0 * j / size;
            //coord_values[1][i * size + j] = -89.0 + 178.0 * i / size;
            coord_values[1][i * size + j] = -89.0 + 178.0 * i / size;
            //coord_values[0][i * size + j] = fRand(0.0, 359.0);
            //coord_values[1][i * size + j] = fRand(-89.0, 89.0);
        } 

    coord_values[0][0] = 0.0;
    coord_values[1][0] = -90.0;
    coord_values[0][299] = 0.0;
    coord_values[1][299] = 90.0;
}


Grid_info_manager::Grid_info_manager(double *coord[2], int num)
{

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
    //*min_lat = -89.0;
    //*max_lat =  89.0;
    *min_lat = -90.0;
    *max_lat =  90.0;
    *min_lon =   0.0;
    //*max_lon = 359.0;
    *max_lon = 360.0;
}

bool Grid_info_manager::is_grid_cyclic(int grid_id)
{
    return true;
}
int Grid_info_manager::get_polar_points(int grid_id, char polar)
{
    return 1;
}
